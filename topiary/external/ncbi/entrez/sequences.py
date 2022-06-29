"""
Use entrez to download protein sequences from the NCBI.
"""

import topiary
from topiary._private import check

from Bio import Entrez

from tqdm.auto import tqdm

import os, urllib, http, time
import multiprocessing as mp

def _get_sequences_thread(args):
    """
    Download ids from the NCBI. args is interpreted as:

    Parameters
    ----------
        args: list of arguments controlling entrez query. Expanded to:
            index: number of request, used for sorting results after multithreading
            ids: string with comma-separated list of ids to download
            num_tries_allowed: how many attempts should be made to actually do
                               a given query
            lock: lock to control number of requests per thread per second
            queue: multiprocessing queue in which to store results
    Return
    ------
        None
    """

    index = args[0]
    ids = args[1]
    num_tries_allowed = args[2]
    lock = args[3]
    queue = args[4]

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

    # Record out and initial index in the output queue.
    queue.put((index,out))

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
    fasta_file : str
        file holding all downloaded sequences.
    """

    block_size = check.check_int(block_size,
                                             "block_size",
                                             minimum_allowed=1)
    num_tries_allowed = check.check_int(num_tries_allowed,
                                                    "num_tries_allowed",
                                                    minimum_allowed=1)
    num_threads = check.check_int(num_threads,
                                              "num_threads",
                                              minimum_allowed=-1)
    if num_threads == 0:
        err = "\nnum_threads cannot be zero. It can be -1 (use all available),\n"
        err += "or any integer > 0.\n\n"
        raise ValueError(err)


    # Figure out number of threads to use
    if num_threads < 0:
        try:
            num_threads = mp.cpu_count()
        except NotImplementedError:
            num_threads = os.cpu_count()
            if num_threads is None:
                print("Could not determine number of cpus. Using single thread.\n",flush=True)
                num_threads = 1

    # Figure out how many blocks we're going to download
    num_blocks = len(to_download) // block_size
    if len(to_download) % block_size > 0:
        num_blocks += 1
    print(f"Downloading {num_blocks} blocks of ~{block_size} sequences... ",flush=True)

    # queue will hold results from each download batch.
    manager = mp.Manager()
    queue = manager.Queue()
    lock = manager.Lock()
    with mp.Pool(num_threads) as pool:

        # This is a bit obscure. Build a list of args to pass to the pool. Each
        # tuple of args matches the args in _get_sequences_thread.
        # all_args has all len(df) recip blast runs we want to do.
        all_args = []
        for i in range(0,len(to_download),block_size):
            ids = ",".join(to_download[i:(i+block_size)])
            all_args.append((i,ids,num_tries_allowed,lock,queue))

        # Black magic. pool.imap() runs a function on elements in iterable,
        # filling threads as each job finishes. (Calls _get_sequences_thread
        # on every args tuple in all_args). tqdm gives us a status bar.
        # By wrapping pool.imap iterator in tqdm, we get a status bar that
        # updates as each thread finishes.
        list(tqdm(pool.imap(_get_sequences_thread,all_args),
                  total=len(all_args)))

    # Get results out of the queue. Sort by i to put in correct order
    results = []
    while not queue.empty():
        results.append(queue.get())
    results.sort()

    # Get final downloaded stuff
    total_download = []
    for r in results:
        total_download.append(r[1])

    print("Done.")

    return "".join(total_download)
