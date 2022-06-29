"""
Run BLAST against a remote NCBI database.
"""

import topiary
from topiary._private import check
from topiary._private import animation
from topiary._private import threads

from .util import read_blast_xml, _standard_blast_args_checker

from Bio import Entrez
from Bio.Blast import NCBIXML, NCBIWWW

import numpy as np
import pandas as pd

from tqdm.auto import tqdm

import sys, urllib, http, copy, os, time, random, string
import multiprocessing as mp

def _prepare_for_blast(sequence,
                       db,
                       taxid,
                       blast_program,
                       hitlist_size,
                       e_value_cutoff,
                       gapcosts,
                       url_base,
                       kwargs):
    """
    Take inputs to ncbi_blast, check arguments, and do initial processing.
    Outputs are a list of validated sequences, a dictionary of keyword arguments
    to past to qblast for each query, and whether or not to return a single df
    or list of df at the end.

    Parameters
    ----------
        sequence: sequence as a string OR list of string sequences
        db: NCBI blast database
        taxid: taxid for limiting blast search (default to no limit)
        blast_program: NCBI blast program to use (blastp, tblastn, etc.)
        hitlist_size: download only the top hitlist_size hits
        e_value_cutoff: only take hits with e_value better than e_value_cutoff
        gapcost: BLAST gapcosts (length 2 tuple of ints)
        url_base: NCBI base url
        kwargs: extra keyword arguments are passed directly to Bio.Blast.NCBIWWW.qblast,
                overriding anything constructed above. You could, for example, pass
                entrez_query="txid9606[ORGN] or txid10090[ORGN]" to limit blast
                results to hits from human or mouse. This would take precedence over
                any taxid specified above.

    Return
    ------
        sequence_list, blast_kwargs, return_singleton
    """

    # Deal with standard input
    out = _standard_blast_args_checker(sequence,
                                       hitlist_size,
                                       e_value_cutoff,
                                       gapcosts)
    sequence_list = out[0]
    hitlist_size = out[1]
    e_value_cutoff = out[2]
    gapcosts = out[3]
    return_singleton = out[4]

    # -------------------------------------------------------------------------
    # Deal with db input

    db_is_bad = True
    if type(db) is str:
        if len(db) > 1:
            db_is_bad = False

    if db_is_bad:
        err = "\ndb should be a string indicating an NCBI database (e.g. 'nr')\n\n"
        raise ValueError(err)

    # -------------------------------------------------------------------------
    # Deal with taxid input

    taxid_out = []
    taxid_is_bad = True

    # Is it int-like? Convert integer input to list of one string. We're
    # stringent on this -- no floats -- because an int cast will silently
    # round down. If someone used a taxid like 960.6 (extra .) it would be
    # interpreted as 960, which is very much the wrong taxid
    if np.issubdtype(type(taxid),np.integer):
        taxid = [f"{taxid}"]

    # This wacky line sees if the taxid is iterable, but not a type *class*.
    # Catches weird edge case where user passes in str or int as a taxid
    if hasattr(taxid,"__iter__") and type(taxid) is not type:

        taxid_is_bad = False

        # If taxid is a single string, convert it to a list of one string
        if type(taxid) is str:
            taxid = [taxid]

        # Go through list of taxids and put in correct format
        for t in taxid:
            if type(t) is str:
                taxid_out.append(f"txid{t}[ORGN]")
            elif np.issubdtype(type(t),np.integer):
                taxid_out.append(f"txid{t}[ORGN]")
            else:
                taxid_is_bad = True
                break

    # If taxid was None to begin with, ignore it
    else:
        if taxid is None:
            taxid_out = []
            taxid_is_bad = False

    if taxid_is_bad:
        err = "\ntaxid should be either a single ncbi taxon id (e.g. 9606 for\n"
        err += "human) or list of ncbi taxon ids. Individual ids can either be\n"
        err += "int or str. These are passed to NCBI without further validation.\n\n"
        raise ValueError(err)

    # -------------------------------------------------------------------------
    # Deal with blast_program

    blast_program_is_bad = True
    if type(blast_program) is str:
        if len(blast_program) > 1:
            blast_program_is_bad = False

    if blast_program_is_bad:
        err = "\nblast_program should be a string identifier for the blast\n"
        err += "program to use (i.e. blastp, tblastn, etc.)\n"
        raise ValueError(err)

    # -------------------------------------------------------------------------
    # Deal with url_base

    url_base_is_bad = True
    if type(url_base) is str:
        if len(url_base) > 1:
            url_base_is_bad = False

    if url_base_is_bad:
        err = "\nurl_base should be a string holding a url pointing to the\n"
        err += "blast server you wish to use.\n"
        raise ValueError(err)


    # -------------------------------------------------------------------------
    # Construct qblast keywords

    # keyword arguments to pass to qblast *besides* sequence
    blast_kwargs = {"program":blast_program,
                    "database":db,
                    "hitlist_size":f"{hitlist_size}",
                    "expect":f"{e_value_cutoff}",
                    "gapcosts":f"{gapcosts[0]} {gapcosts[1]}",
                    "url_base":url_base}

    # Construct taxid entrez_query
    if len(taxid_out) > 0:
        blast_kwargs["entrez_query"] = " or ".join(taxid_out)

    # Capture the rest of kwargs, overwriting anything automatically made
    for k in kwargs:
        blast_kwargs[k] = kwargs[k]

    return sequence_list, blast_kwargs, return_singleton

def _construct_args(sequence_list,
                    blast_kwargs,
                    max_query_length,
                    num_tries_allowed,
                    num_threads=-1,
                    test_num_cores=None):
    """
    Construct a list of arguments to pass to each thread that will run a
    blast query.

    Parameters
    ----------
    sequence_list : list
        list of sequences as strings
    blast_kwargs : dict
        keyword arguments to pass to blast call
    max_query_length : int
        maximum string length accepted by the server. if the query is too long,
        this function will break it into multiple requests, each sent to ncbi.
    num_tries_allowed : int
        try num_tries_allowed times in case of timeout
    num_threads : int, default=-1
        number of threads to use (locally). if -1 use all available.
    test_num_cores : int
        send in a hacked number of cores for testing purposes

    Returns
    ------
    kwargs_list : list
        list of kwargs to pass for each calculation
    num_threds : int
        number of threads to use
    """

    # Validate inputs that have not yet been validated.
    max_query_length = check.check_int(max_query_length,
                                       "max_query_length",
                                       minimum_allowed=1)

    num_tries_allowed = check.check_int(num_tries_allowed,
                                        "num_tries_allowed",
                                        minimum_allowed=1)

    num_threads = threads.get_num_threads(num_threads,test_num_cores)

    # Make sure sequence_list is an array
    sequence_list = np.array(sequence_list)

    # Break sequences up into blocks
    num_sequences = len(sequence_list)
    block_size = num_sequences//num_threads

    # If more threads than sequences, block_size is 1. Only run one per thread.
    if block_size < 1:
        block_size = 1
        num_threads = num_sequences

    # Windows will hold indexes for each block
    windows = [block_size for _ in range(num_sequences//block_size)]
    remainder = num_sequences % block_size
    if remainder/block_size > 0.5:
        windows.append(remainder)
    else:
        counter = 0
        while remainder > 0:
            windows[counter] += 1
            remainder -= 1
    windows.insert(0,0)

    windows = np.cumsum(windows)

    # Blocks will allow us to tile over all sequences
    split_sequences = []
    counter = 0
    for i in range(len(windows)-1):
        split_sequences.append([])
        for seq in sequence_list[windows[i]:windows[i+1]]:

            new_sequence = f">count{counter}\n{seq}\n"

            # If you have a super, super long sequence...
            if len(new_sequence) > max_query_length:
                err = "\nat least one sequence is, by itself, greater than\n"
                err += f"the maximum server query length '{max_query_length}'.\n"
                err += "To fix, either trim this sequence, remove it from the\n"
                err += "query list, increase the max_query_length (not\n"
                err += "recommended), or use a local blast database rather\n"
                err += "than the NCBI database. Sequence is:\n\n"
                err += f"{seq}\n\n"
                raise ValueError(err)

            new_length = len(new_sequence)
            current_length = len("".join(split_sequences[-1]))
            if current_length + new_length > max_query_length:
                split_sequences.append([new_sequence])
            else:
                split_sequences[-1].append(new_sequence)

            counter += 1

    final_sequences = []
    for s in split_sequences:
        final_sequences.append("".join(s))

    kwargs_list = []
    for i, seq in enumerate(final_sequences):
        query = copy.deepcopy(blast_kwargs)
        query["sequence"] = seq
        kwargs_list.append({"this_query":query,
                            "num_tries_allowed":num_tries_allowed})

    return kwargs_list, num_threads


def _ncbi_blast_thread_function(this_query,num_tries_allowed,lock):
    """
    Run an NCBIWWW.qblast call on a single thread, making several attempts.
    Put results in a multiprocessing queue.

    Parameters
    ----------
    this_query : dict
        kwargs to pass to NCBIWWW.qblast
    num_tries : int
        number of tries before timing out
    lock : multiprocessing.Manager().Lock()
        lock used to prevent hammering NCBI servers at too high of a rate.

    Returns
    ------
    out_df : pandas.DataFrame
        dataframe with BLAST hits
    """

    # While we haven't run out of tries...
    tries = 0
    while tries < num_tries_allowed:

        # Try to read output.
        try:

            # NCBI limits requests to 3 per second for normal users. Use a lock
            # when launching that sleeps for 0.5 seconds. This means we'll make
            # a maximum of two requests per second across all threads and
            # should avoid the NCBI limit.
            lock.acquire()
            try:
                time.sleep(0.5)
            finally:
                lock.release()

            result = NCBIWWW.qblast(**this_query)

            # Write output to an xml file. NCBI can spit out trashed XML
            # with CREATE_VIEW\n\n\n randomly injected between <hit> entries.
            # biopython chokes on the input. To fix this, pull down the XML,
            # clean up, write to a file, then pass the file handle back to
            # biopython.

            # Get rid of nastiness if present.
            contents = result.readlines()
            contents = [c for c in contents if c.strip() not in ["","CREATE_VIEW"]]

            # Write temporary file
            tmp_root = "".join([random.choice(string.ascii_letters)
                                for _ in range(10)])
            tmp_file = f"{tmp_root}_ncbi-blast-result.xml"
            f = open(tmp_file,"w")
            f.write("".join(contents))
            f.close()

            # Read output xml file and parse
            f = open(tmp_file,"r")

            # Clean up
            p = NCBIXML.parse(f)
            out = []
            for r in p:
                out.append(r)
            f.close()

            # If parsing successful, nuke temporary file
            os.remove(tmp_file)

        # If some kind of http error or timeout, set out to None
        except (urllib.error.URLError,urllib.error.HTTPError,http.client.IncompleteRead):
            out = None

        # If out is None, try again. If not, break out of loop--success!
        if out is None:
            tries += 1
            continue
        else:
            break

    # We didn't get result even after num_tries_allowed tries. Throw
    # an error.
    if out is None:
        err = "\nProblem accessing with NCBI server. We got no output after\n"
        err += f"{num_tries_allowed} attempts. This is likely due to a server\n"
        err += "timeout. Try again at a later time.\n"
        raise RuntimeError(err)

    # Parse output
    out_df = read_blast_xml(out)

    return out_df


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

    all_df = []
    for hit in hits:

        # Break big dataframe into a list of dataframes, one for each query sequence
        queries = np.unique(hit["query"])

        # Trim off "count" from "countNUMBER" queries and use to sort
        query_order = [(q[5:],q) for q in queries]
        query_order.sort()

        # Go through all queries in this hit dataframe
        for q in query_order:

            # Get dataframe
            this_df = hit.loc[hit["query"] == q[1],:]

            # No hits, return empty dataframe
            if len(this_df) == 1 and pd.isna(this_df["accession"].iloc[0]):
                all_df.append(pd.DataFrame())

            # Record hits
            else:
                all_df.append(this_df)


    # Singleton hit -- pop out of list
    if return_singleton:
        all_df = all_df[0]

    return all_df


def ncbi_blast(sequence,
               db="nr",
               blast_program="blastp",
               taxid=None,
               hitlist_size=100,
               e_value_cutoff=0.001,
               gapcosts=(11,1),
               url_base="https://blast.ncbi.nlm.nih.gov/Blast.cgi",
               max_query_length=80000,
               num_tries_allowed=5,
               num_threads=1,
               verbose=False,
               **kwargs):
    """
    Perform a blast query against a remote NCBI blast database. Takes a sequence
    or list of sequences and returns a list of topiary dataframes containing hits
    for each sequence.

    Parameters
    ----------
    sequence : str or list
        sequence as a string OR list of string sequences
    db : str, default="nr"
        name of ncbi blast database
    blast_program : str, default="blastp"
        NCBI blast program to use (i.e. "blastp", "tblastn", etc.)
    taxid : str or int or list or None, default=None
        limit ncbi blast search to specified taxid. If None, do not limit
        search. If single str/int, interpret as an integer taxid (i.e. 9606 for
        Homo sapiens). If list of str/int, interpret as multiple taxid (i.e.
        [9606,10090] for Homo sapiens OR Mus musculus).
    hitlist_size : int, default=100
        return only the top hitlist_size hits
    e_value_cutoff : float, default=0.001
        only return hits with e_value better than e_value_cutoff
    gapcosts : tuple, default=(11,1)
        BLAST gapcosts (length 2 tuple of ints)
    url_base : str, default="https://blast.ncbi.nlm.nih.gov/Blast.cgi"
        NCBI base url
    max_query_length : int, default=80000
        maximum string length accepted by the server. if the query is longer
        than max_query_length, the function will break it into multiple requests
        shorter than max_query_length, then combine the results.
    num_tries_allowed : int, default=5
        try num_tries_allowed times in case of server timeout
    num_threads : int, default=1
        number of threads to use. if -1, use all available.
    verbose : bool, default=False
        whether or not to use verbose output
    **kwargs : dict, optional
        extra keyword arguments are passed directly to Bio.Blast.NCBIWWW.qblast,
        overriding anything constructed above. You could, for example, pass
        entrez_query="txid9606[ORGN] or txid10090[ORGN]" to limit blast
        results to hits from human or mouse. This would take precedence over
        any taxid specified above.

    Returns
    -------
    blast_output : pandas.DataFrame or list
        If a single sequence is passed in, return a dataframe of hits. If a list
        of queries is passed in, returns a list of dataframes--one entry per
        query. If no hits are found for a given query, return an empty dataframe.
    """

    # Generate a list of sequences, keywords to pass to qblast, and whether or
    # not to return a single or list of dataframes
    prep = _prepare_for_blast(sequence=sequence,
                              db=db,
                              taxid=taxid,
                              blast_program=blast_program,
                              hitlist_size=hitlist_size,
                              e_value_cutoff=e_value_cutoff,
                              gapcosts=gapcosts,
                              url_base=url_base,
                              kwargs=kwargs)

    sequence_list = prep[0]
    blast_kwargs = prep[1]
    return_singleton = prep[2]

    # Get arguments to pass to each thread
    kwargs_list, num_threads = _construct_args(sequence_list=sequence_list,
                                               blast_kwargs=blast_kwargs,
                                               max_query_length=max_query_length,
                                               num_threads=num_threads,
                                               num_tries_allowed=num_tries_allowed)

    # Run multi-threaded blast
    print(f"Performing {len(sequence_list)} BLAST queries against the NCBI {db} database")
    print(f"on {num_threads} threads. Depending on the server load, this could")
    print(f"take awhile. This is a good time to grab a cup of coffee.",flush=True)

    a = animation.WaitingAnimation()
    a.start()

    hits = threads.thread_manager(kwargs_list,
                                  _ncbi_blast_thread_function,
                                  num_threads,
                                  progress_bar=False,
                                  pass_lock=True)

    a.stop()
    print("BLAST query complete.",flush=True)


    # Combine hits into dataframes, one for each query
    out_df = _combine_hits(hits,return_singleton)

    return out_df
