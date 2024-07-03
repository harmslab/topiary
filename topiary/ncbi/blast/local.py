"""
Run BLAST against a local database.
"""

import topiary
from topiary._private import check
from topiary._private import threads
from topiary._private import interface
from .util import _standard_blast_args_checker
from .read import read_blast_xml

#import Bio.Blast.Applications as apps

import numpy as np
import pandas as pd

import os, string, random, subprocess, copy

def _prepare_for_blast(sequence,
                       db,
                       blast_program,
                       hitlist_size,
                       e_value_cutoff,
                       gapcosts,
                       kwargs,
                       test_skip_blast_program_check=False):
    """
    Take inputs to local_blast, check arguments, and do initial processing.
    Outputs are a list of validated sequences, the blast function to use,
    a dictionary of keyword arguments to past to blast for each query, and
    whether or not to return a single df or list of df at the end.

    Parameters
    ----------
    sequence : str or list
        sequence as a string OR list of string sequences
    db: str
        name of local blast database
    blast_program : str
        NCBI blast program to use (blastp, tblastn, etc.)
    hitlist_size : int
        download only the top hitlist_size hits
    e_value_cutoff : float
        only take hits with e_value better than e_value_cutoff
    gapcosts : tuple
        BLAST gapcosts (length 2 tuple of ints)
    kwargs : dict
        extra keyword arguments are passed directly to blast. keys should be 
        valid arguments to pass to blast program {"db":"blah"}. This function
        will append '-' to front of those arguments. 
    test_skip_blast_program_check : bool
        skip the check for a working blast program (for testing)

    Return
    ------
    sequence_list : list
        list of sequences
     blast_args : list
        list containing the blast command to run (e.g. ["blastp","-db","blah"...]),
    return_singleton : bool
        whether or not to return final output as a single sequence or list of
        sequences.
    """

    recognized_functions = ["blastp",
                            "blastn",
                            "blastx",
                            "tblastn",
                            "tblastx",
                            "psiblast",
                            "rpsblast",
                            "rpstblastn",
                            "deltablast"]

    if not os.path.exists(f"{db}.psq"):
        err = f"db {db}.psq not found!\n"
        raise FileNotFoundError(err)

    # Make sure we can actually run the local blasting
    if blast_program not in recognized_functions:
        err = "\nblast_program '{}' not recognized.\n\nAllowed programs:\n".format(blast_program)
        for f in recognized_functions:
            err += "    {}\n".format(f)
        raise ValueError(err)

    # Make sure the blast program is installed.
    if not test_skip_blast_program_check:
        try:
            _ = subprocess.run([blast_program],stderr=subprocess.DEVNULL)
        except FileNotFoundError:
            err = f"\nBLAST program {blast_program} is not in path. Is it installed?\n\n"
            raise FileNotFoundError(err)

    out = _standard_blast_args_checker(sequence,
                                       hitlist_size,
                                       e_value_cutoff,
                                       gapcosts)
    sequence_list = out[0]
    hitlist_size = out[1]
    e_value_cutoff = out[2]
    gapcosts = out[3]
    return_singleton = out[4]

    # Construct blast arguments to pass in as a process
    blast_args = [blast_program]
    blast_args.extend(["-db",db])
    blast_args.extend(["-outfmt","5"])
    blast_args.extend(["-max_target_seqs",f"{hitlist_size:d}"])
    blast_args.extend(["-threshold",f"{e_value_cutoff:.5e}"])
    blast_args.extend(["-gapopen",f"{gapcosts[0]:d}"])
    blast_args.extend(["-gapextend",f"{gapcosts[1]:d}"])
    
    # Load in any keyword arguments
    for k in kwargs:
        blast_args.extend([f"-{k}","{}".format(kwargs[k])])

    return sequence_list, blast_args, return_singleton

def _construct_args(sequence_list,
                    blast_args,
                    keep_blast_xml=False,
                    block_size=20,
                    num_threads=-1,
                    manual_num_cores=None):
    """
    Construct a list of arguments to pass to each thread in the pool.

    Parameters
    ----------
    sequence_list : list
        list of sequences as strings
    blast_args : list
        list of blast arguments
    keep_blast_xml : bool, default=False
        whether or not to keep temporary files
    block_size : int, default=20
        break into block_size sequence chunks
    num_threads : int, default=-1
        number of threads to use. if -1, use all available.
    manual_num_cores : int, optional
        use exactly this number of cores
        
    Returns
    -------
    kwargs_dict : dict
        dictionary of keyword arguments to pass
    num_threads : int
        number of threads to use for the calculation
    """

    # Validate inputs that have not yet been validated.
    block_size = check.check_int(block_size,
                                 "block_size",
                                 minimum_allowed=1)

    num_threads = threads.get_num_threads(num_threads,manual_num_cores)

    keep_blast_xml = check.check_bool(keep_blast_xml,"keep_blast_xml")

    # Determine number of threads useful for this problem. It's not worth
    # chopping up a super small set of comparisons
    max_useful_threads = len(sequence_list)//block_size
    if max_useful_threads < 1:
        max_useful_threads = 1
    if max_useful_threads > num_threads:
        max_useful_threads = num_threads

    # Set number of threads
    if num_threads > max_useful_threads:
        num_threads = max_useful_threads

    # Break sequences up into blocks
    num_sequences = len(sequence_list)

    num_windows = num_sequences//block_size
    if num_windows == 0:
        num_windows = 1
        block_size = num_sequences

    windows = [block_size for _ in range(num_windows)]
    remainder = num_sequences % block_size

    if remainder/block_size >= 0.5:
        windows.append(remainder)
    else:
        counter = 0
        while remainder > 0:
            windows[counter] += 1
            remainder -= 1
            counter += 1
            if counter >= len(windows):
                counter = 0

    windows.insert(0,0)
    windows = np.cumsum(windows)

    # Blocks will allow us to tile over whole sequence
    kwargs_list = []
    for i in range(len(windows)-1):

        i_block = (windows[i],windows[i+1])

        kwargs_list.append({"sequence_list":sequence_list,
                            "index":i_block,
                            "blast_args":blast_args,
                            "keep_blast_xml":keep_blast_xml})


    return kwargs_list, num_threads


def _local_blast_thread_function(sequence_list,
                                 index,
                                 blast_args,
                                 keep_blast_xml):
    """
    Run local blast on a list of sequences.

    Parameters
    ----------
    sequence_list : list
        list of sequences
    index : tuple
        indexes to pull from sequence_list
    blast_args : list
        list holding blast command to past (e.g. ["blastp","-db","yo",...])
    keep_blast_xml : bool
        whether or not to keep temporary files

    Returns
    -------
    out_df : pandas.DataFrame
        dataframe containing blast hits
    """

    # make a 10-character random string for temporary files
    tmp_file_root = "".join([random.choice(string.ascii_letters) for i in range(10)])
    input_file = "topiary-tmp_{}_blast-in.fasta".format(tmp_file_root)
    out_file = "topiary-tmp_{}_blast-out.xml".format(tmp_file_root)

    f = open(input_file,'w')
    for i in range(index[0],index[1]):
        f.write("".join(f">count{i}\n{sequence_list[i]}\n"))
    f.close()

    blast_args = copy.deepcopy(blast_args)
    blast_args.extend(["-query",input_file])
    blast_args.extend(["-out",out_file])

    interface.launch(blast_args,
                     run_directory=".",
                     suppress_output=True)

    


    # Parse output
    try:
        out_dfs, _ = read_blast_xml(out_file)
    except FileNotFoundError:
        err = "\nLocal blast failed on sequence:\n"
        err += f"    '{sequence_list[i]}'\n\n"
        raise RuntimeError(err)

    # If parsing successful, nuke temporary file
    if not keep_blast_xml:
        os.remove(input_file)
        os.remove(out_file)

    return out_dfs[0]



def _combine_hits(hits,return_singleton):
    """
    Parse a list of hits resulting from a set of blast queries and return a
    list of dataframes, one for each query. Results will be ordered first by the
    order they occur in the list, then sorted by query. (Assumes queries have
    the form countNUMBER when it does sorting.

    Parameters
    ----------
        hits: list of dataframes containing hits
        return_singleton: bool. whether or not to return a single dataframe or
                          list of dataframes

    Return
    ------
        combined blast hits (list of dataframes or single dataframe)
    """

    # Construct a list of output dataframes, one for each query sequence
    out_df = []
    for h in hits:

        queries = np.unique(h["query"])
        queries = [(int(q[5:]),q) for q in queries]
        queries.sort()
        queries = [q[1] for q in queries]

        for q in queries:
            this_df = h.loc[h["query"] == q,:]

            # No hits, return empty dataframe
            if len(this_df) == 1 and pd.isna(this_df["accession"].iloc[0]):
                out_df.append(pd.DataFrame())

            # Record hits
            else:
                out_df.append(this_df)

    # Singleton hit -- pop out of list
    if return_singleton:
        out_df = out_df[0]

    return out_df


def local_blast(sequence,
                db,
                blast_program="blastp",
                hitlist_size=100,
                e_value_cutoff=0.001,
                gapcosts=(11,1),
                keep_blast_xml=False,
                num_threads=-1,
                block_size=20,
                **kwargs):
    """
    Perform a blast query against a local blast database. Takes a sequence or
    list of sequences and returns a list of topiary dataframes containing hits
    for each sequence.

    Parameters
    ----------
    sequence : str or list
        sequence as a string OR list of string sequences
    db : str
        name of local blast database
    blast_program : str, default="blastp"
        NCBI blast program to use (i.e. "blastp", "tblastn", etc.)
    hitlist_size : int, default=100
        return only the top hitlist_size hits
    e_value_cutoff : float, default=0.001
        only return hits with e_value better than e_value_cutoff
    gapcosts : tuple, default=(11,1)
        BLAST gapcosts (length 2 tuple of ints)
    keep_blast_xml : bool, default=False
        whether or not to keep temporary blast xml files
    num_threads : int, default=-1
        number of threads to use. if -1, use all available.
    block_size : int, default=20
        run blast in blocks of block_size sequences
    **kwargs : dict, optional
        extra keyword arguments are passed directly to the blast call. keys
        should be valid arguments to pass to a blast program. For example, 
        this could be {"db":"blah"}. The function will append "-" to the front,
        making "-db blah" when calling the blast program. 
        
    Returns
    -------
    blast_output : pandas.DataFrame or list
        If a single sequence is passed in, return a dataframe of hits. If a list
        of queries is passed in, returns a list of dataframes--one entry per
        query. If no hits are found for a given query, return an empty dataframe.
    """

    # Generate BLAST input
    prep = _prepare_for_blast(sequence=sequence,
                              db=db,
                              blast_program=blast_program,
                              hitlist_size=hitlist_size,
                              e_value_cutoff=e_value_cutoff,
                              gapcosts=gapcosts,
                              kwargs=kwargs)

    sequence_list = prep[0]
    blast_args = prep[1]
    return_singleton = prep[2]

    # Construct a list of arguments to pass into _thread_manager
    kwargs_list, num_threads = _construct_args(sequence_list=sequence_list,
                                               blast_args=blast_args,
                                               keep_blast_xml=keep_blast_xml,
                                               block_size=block_size,
                                               num_threads=num_threads)

    # Run multi-threaded local blast
    hits = threads.thread_manager(kwargs_list,
                                  _local_blast_thread_function,
                                  num_threads)

    # Combine hits into dataframes, one for each query
    out_df = _combine_hits(hits,return_singleton)

    return out_df
