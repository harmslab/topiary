__description__ = \
"""
Run BLAST against a local database.
"""
__author__ = "Michael J. Harms"
__date__ = "2021-04-08"

import topiary
from topiary import _arg_processors
from .util import read_blast_xml, _standard_blast_args_checker

import Bio.Blast.Applications as apps

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

import sys, os, string, random, subprocess, copy
import multiprocessing as mp

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
        sequence: sequence as a string OR list of string sequences
        db: name of local blast database
        blast_program: NCBI blast program to use (blastp, tblastn, etc.)
        hitlist_size: download only the top hitlist_size hits
        e_value_cutoff: only take hits with e_value better than e_value_cutoff
        gapcosts: BLAST gapcosts (length 2 tuple of ints)
        kwargs: extra keyword arguments are passed directly to
                apps.NcbiblastXXXCommandline, overriding anything constructed
                above.
        test_skip_blast_program_check: skip the check for a working blast
                                       program (for testing)

    Return
    ------
        dataframe (single query) or list of dataframes (multiple sequences).
        if no hits found, warns and returns None.
    """

    recognized_functions = {"blastp":apps.NcbiblastpCommandline,
                            "blastn":apps.NcbiblastnCommandline,
                            "blastx":apps.NcbiblastxCommandline,
                            "tblastn":apps.NcbitblastnCommandline,
                            "tblastx":apps.NcbitblastxCommandline,
                            "psiblast":apps.NcbipsiblastCommandline,
                            "rpsblast":apps.NcbirpsblastCommandline,
                            "rpstblastn":apps.NcbirpstblastnCommandline,
                            "deltablast":apps.NcbideltablastCommandline}

    if not os.path.exists(f"{db}.psq"):
        err = f"db {db}.psq not found!\n"
        raise FileNotFoundError(err)

    # Make sure we can actually run the local blasting
    try:
        blast_function = recognized_functions[blast_program]
    except (KeyError,TypeError):
        err = "\nblast_program '{}' not recognized.\n\nAllowed programs:\n".format(blast_program)
        for k in recognized_functions.keys():
            err += "    {}\n".format(k)
        raise ValueError(err)

    # Make sure the blast program is installed.
    if not test_skip_blast_program_check:
        try:
            ret = subprocess.run([blast_program],stderr=subprocess.DEVNULL)
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

    gaps = '{} {}'.format(*gapcosts)

    # Construct keyword arguments to pass to function
    blast_kwargs = {"cmd":blast_program,
                    "db":db,
                    "outfmt":5,
                    "max_target_seqs":hitlist_size,
                    "threshold":e_value_cutoff,
                    "gapopen":gapcosts[0],
                    "gapextend":gapcosts[1]}
    for k in kwargs:
        blast_kwargs[k] = kwargs[k]


    return sequence_list, blast_function, blast_kwargs, return_singleton

def _construct_args(sequence_list,
                    blast_function,
                    blast_kwargs,
                    keep_tmp=False,
                    block_size=20,
                    num_threads=-1,
                    test_num_cores=None):
    """
    Construct a list of arguments to pass to each thread in the pool.

    Parameters
    ----------
        sequence_list: list of sequences as strings
        blast_function: blast function to use
        blast_kwargs: keyword arguments to pass to blast call
        keep_tmp: whether or not to keep temporary files
        num_threads: number of threads to use. if -1, use all available.
        block_size: break into block_size sequence chunks

    Return
    ------
        list of args to pass for each calculation, number of threads
    """

    # Validate inputs that have not yet been validated.
    block_size = _arg_processors.process_int(block_size,
                                             "block_size",
                                              minimum_allowed=1)

    num_threads = _arg_processors.process_int(num_threads,
                                              "num_threads",
                                              minimum_allowed=-1)
    if num_threads == 0:
        err = "\nnum_threads cannot be zero. It can be -1 (use all available),\n"
        err += "or any integer > 0.\n\n"
        raise ValueError(err)

    keep_tmp = _arg_processors.process_bool(keep_tmp,"keep_tmp")

    # Try to figure out how many cores are available
    if test_num_cores is not None:
        num_cores = test_num_cores
    else:
        try:
            num_cores = mp.cpu_count()
        except NotImplementedError:
            num_cores = os.cpu_count()

        # If we can't figure it out, revert to 1
        if num_cores is None:
            print("Could not determine number of cpus. Using single thread.",flush=True)
            num_cores = 1

    # Determine number of threads useful for this problem. It's not worth
    # chopping up a super small set of comparisons
    max_useful_threads = len(sequence_list)//block_size
    if max_useful_threads < 1:
        max_useful_threads = 1
    if max_useful_threads > num_cores:
        max_useful_threads = num_cores

    # Set number of threads
    if num_threads == -1 or num_threads > max_useful_threads:
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

    windows.insert(0,0)
    windows = np.cumsum(windows)

    # Blocks will allow us to tile over whole sequence
    all_args = []
    for i in range(len(windows)-1):

        i_block = (windows[i],windows[i+1])

        all_args.append([sequence_list,
                         i_block,
                         blast_function,
                         blast_kwargs,
                         keep_tmp])

    return all_args, num_threads


def _thread_manager(all_args,num_threads):
    """
    Run a bunch of blast jobs in a mulithreaded fashion. Should only be called
    by local_blast.

    Parameters
    ----------
        all_args: list of args to pass for each calculatio
        nm_threads: number of threads to use

    Return
    ------
        list of dataframes with blast results
    """

    print(f"Performing {len(all_args)} blast queries on {num_threads} threads.",
          flush=True)

    # queue will hold results from each run. Append to each entry in all_args
    queue = mp.Manager().Queue()
    for i in range(len(all_args)):
        all_args[i].append(queue)

    with mp.Pool(num_threads) as pool:

        # Black magic. pool.imap() runs a function on elements in iterable,
        # filling threads as each job finishes. (Calls _blast_thread
        # on every args tuple in all_args). tqdm gives us a status bar.
        # By wrapping pool.imap iterator in tqdm, we get a status bar that
        # updates as each thread finishes.
        list(tqdm(pool.imap(_thread,all_args),total=len(all_args)))

    # Get results out of the queue.
    results = []
    while not queue.empty():
        results.append(queue.get())

    # Sort results
    results.sort()
    hits = [r[1] for r in results]

    return hits

def _thread(args):
    """
    Run blast on a thread. Should only be called via _thread_manager.
    Puts resulting hits as a pandas dataframe into the queue

    Parameters
    ----------
        args. list that is expanded as follows:

        sequence_list: list of sequences as strings
        index: sequence to grab
        blast_kwargs: keyword arguments to pass to blast call
        blast_function: blast function to call
        keep_tmp: whether or not to keep temporary files
        queue: multiprocessing queue for storing results

    Return
    ------
        None
    """

    # parse args
    sequence_list = args[0]
    index = args[1]
    blast_function = args[2]
    blast_kwargs = args[3]
    keep_tmp = args[4]
    queue = args[5]

    # make a 10-character random string for temporary files
    tmp_file_root = "".join([random.choice(string.ascii_letters) for i in range(10)])
    input_file = "topiary-tmp_{}_blast-in.fasta".format(tmp_file_root)
    out_file = "topiary-tmp_{}_blast-out.xml".format(tmp_file_root)

    f = open(input_file,'w')
    for i in range(index[0],index[1]):
        f.write("".join(f">count{i}\n{sequence_list[i]}\n"))
    f.close()

    blast_kwargs = copy.deepcopy(blast_kwargs)
    blast_kwargs["query"] = input_file
    blast_kwargs["out"] = out_file

    blast_function(**blast_kwargs)()

    # Parse output
    try:
        out_df = read_blast_xml(out_file)
    except FileNotFoundError:
        err = "\nLocal blast failed on sequence:\n"
        err += f"    '{sequence_list[i]}'\n\n"
        raise RuntimeError(err)

    # Delete temporary files
    if not keep_tmp:
        os.remove(input_file)
        os.remove(out_file)

    # update queue
    queue.put((index[0],out_df))

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
                e_value_cutoff=0.01,
                gapcosts=(11,1),
                keep_tmp=False,
                num_threads=-1,
                block_size=20,
                **kwargs):
    """
    Perform a blast query against a local blast database. Takes a sequence or
    list of sequences and returns a list of topiary dataframes containing hits
    for each sequence.

    Parameters
    ----------
        sequence: sequence as a string OR list of string sequences
        db: name of local blast database
        blast_program: NCBI blast program to use (blastp, tblastn, etc.)
        hitlist_size: download only the top hitlist_size hits
        e_value_cutoff: only take hits with e_value better than e_value_cutoff
        gapcosts: BLAST gapcosts (length 2 tuple of ints)
        keep_tmp: whether or not to keep temporary blast output
        num_threads: number of threads to use. if -1, use all available.
        block_size: run blast in blocks of block_size sequences
        kwargs: extra keyword arguments are passed directly to apps.NcbiblastXXXCommandline,
                overriding anything constructed above.

    Return
    ------
        dataframe (single query) or list of dataframes (multiple sequences).
        if no hits found, warns and returns None.
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
    blast_function = prep[1]
    blast_kwargs = prep[2]
    return_singleton = prep[3]

    # Construct a list of arguments to pass into _thread_manager
    all_args, num_threads = _construct_args(sequence_list=sequence_list,
                                            blast_function=blast_function,
                                            blast_kwargs=blast_kwargs,
                                            keep_tmp=keep_tmp,
                                            block_size=block_size,
                                            num_threads=num_threads)

    # Run multi-threaded local blast
    hits = _thread_manager(all_args,num_threads)

    # Combine hits into dataframes, one for each query
    out_df = _combine_hits(hits,return_singleton)

    return out_df
