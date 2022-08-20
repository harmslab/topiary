"""
Use entrez to download protein sequences from the NCBI.
"""

import topiary
from topiary._private import check, threads

from Bio import Entrez, SeqIO

import urllib, http, time, io

def _get_sequences_thread_function(ids,num_tries_allowed,lock):
    """
    Download ids from the NCBI.

    Parameters
    ----------
    ids : str
        string with comma-separated list of ids to download
    num_tries_allowed : int
        how many attempts should be made to actually do a given query
    lock : multiprocessing.Manager().Lock()
        lock used to prevent hammering NCBI servers at too high of a rate.

    Returns
    -------
    out : handle.read()
        output from Bio.Entrez
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

            handle = Entrez.efetch(db="protein",
                                   id=ids,
                                   rettype="fasta",
                                   retmode="text")
            out = handle.read()
            handle.close()

        # If some kind of http error or timeout, set out to None
        except (urllib.error.HTTPError,http.client.IncompleteRead):
            out = None

        # If out is None, try again. If not, break out of loop--success!
        if out is None:
            tries += 1
            continue
        else:
            break

    # We failed to download. Throw an error.
    if out is None:
        err = "could not download sequences\n"
        raise ValueError(err)

    # Extract sequences from the downloaded data
    seq_output = []
    for record in SeqIO.parse(io.StringIO(out.strip()), "fasta"):
        seq_id = str(record.id)
        sequence = str(record.seq)
        seq_output.append((seq_id,sequence))

    # Make sure we brought down as many sequences as input ids. If we didn't,
    # figure out which one(s) we are missing and append None for the missing
    # sequences.
    id_list = ids.split(",")
    if len(id_list) != len(seq_output):

        # Create dictionary keying ids we pulled down to their seq_id,sequence
        # entry. We have to slice out the sequence id part to actually do the
        # matching.  (pdb|1STN|A --> 1STN_A, XX|blah.1 --> blah)
        seq_dict = {}
        for s in seq_output:
            seq_id = s[0]
            split_list = seq_id.split(".")[0].split("|")
            if len(split_list) == 1:
                this_id = split_list[0]
            elif len(split_list) == 2:
                this_id = split_list[1]
            elif len(split_list) == 3:
                this_id = "{}_{}".format(*split_list[1:])
            else:
                err = f"\ncould not parse ncbi sequence identifier '{s}'\n\n"
                raise ValueError(err)

            seq_dict[this_id] = (s[0],s[1])

        # Now try to match input to the output. If this does not work, append
        # None.
        final_output = []
        for id in id_list:
            try:
                final_output.append(seq_dict[id])
            except KeyError:
                print(f"No sequence pulled down for query '{id}'",flush=True)
                final_output.append((id,None))

        seq_output = final_output

    return seq_output

def get_sequences(to_download,
                  block_size=50,
                  num_tries_allowed=10,
                  num_threads=-1):
    """
    Use entrez to download protein sequences from the NCBI.

    Parameters
    ----------
    to_download : list
        list of ncbi ids to download
    block_size : int, default=50
        download in chunks this size
    num_tries_allowed : int, default=10
        number of times to try before giving up and throwing an error.
    num_threads : int, default=-1
        number of threads to use. if -1, use all available.

    Returns
    -------
    seq_output : list
        list of tuples of strings. Each tuple looks like (seq_id,sequence)
    """

    # Check to_download. if list is empty, return empty list
    to_download = check.check_iter(to_download,"to_download")
    if len(to_download) == 0:
        return []

    # Check/parse other arguments
    block_size = check.check_int(block_size,
                                 "block_size",
                                 minimum_allowed=1)
    num_tries_allowed = check.check_int(num_tries_allowed,
                                        "num_tries_allowed",
                                        minimum_allowed=1)

    num_threads = threads.get_num_threads(num_threads)

    # Figure out how many blocks we're going to download
    num_blocks = len(to_download) // block_size
    if len(to_download) % block_size > 0:
        num_blocks += 1
    print(f"Downloading {num_blocks} blocks of ~{block_size} sequences... ",flush=True)

    # Construct kwargs to pass to _get_sequences_thread_function
    kwargs_list = []
    for i in range(0,len(to_download),block_size):
        ids = ",".join(to_download[i:(i+block_size)])
        kwargs_list.append({"ids":ids,
                            "num_tries_allowed":num_tries_allowed})

    # Download sequences in multi-threaded fashion
    results = threads.thread_manager(kwargs_list,
                                     _get_sequences_thread_function,
                                     num_threads,
                                     progress_bar=True,
                                     pass_lock=True)

    # Get final downloaded stuff
    seq_output = []
    for r in results:
        seq_output.extend(r)

    return seq_output
