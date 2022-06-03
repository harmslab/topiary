__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
Functions for downloading sequences off entrez in a multithreaded fashion.
"""

from Bio import Entrez
Entrez.email = "DUMMY_EMAIL@DUMMY_URL.COM"

from tqdm.auto import tqdm

import os, urllib, http, time
import multiprocessing as mp

def _entrez_download_thread(args):
    """
    Download ids from the NCBI. args is interpreted as:

    index: number of request, used for sorting results after multithreading
    ids: string with comma-separated list of ids to download
    num_tries_allowed: how many attempts should be made to actually do
                       a given query
    wait_time_between_requests: how many seconds to wait between failed requests.
    queue: multiprocessing queue in which to store results
    """

    index = args[0]
    ids = args[1]
    num_tries_allowed = args[2]
    wait_time_between_requests = args[3]
    queue = args[4]

    # While we haven't run out of tries...
    tries = 0
    while tries < num_tries_allowed:

        # Try to read output.
        try:
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
            time.sleep(wait_time_between_requests)
            continue
        else:
            break

    # We failed to download. Throw an error.
    if out is None:
        err = "could not download sequences\n"
        raise ValueError(err)

    # Record out and initial index in the output queue.
    queue.put((index,out))

def entrez_download(to_download,block_size=50,num_tries_allowed=10,wait_time_between_requests=1,num_threads=-1):
    """
    Download sequences off of entrez, catching errors. This is done in a
    multi-threaded fashion, allowing multiple requests to NCBI at the same
    time.

    to_download: list of ids to download
    block_size: download in chunks this size
    num_tries_allowed: number of times to try before giving up and throwing
                       an error.
    wait_time_between_requests: how many seconds to wait between failed requests.
    num_threads: number of threads to use. if -1, use all available.
    """

    # Figure out number of threads to use
    if num_threads < 0:
        try:
            num_threads = mp.cpu_count()
        except NotImplementedError:
            num_threads = os.cpu_count()
            if num_threads is None:
                print("Could not determine number of cpus. Using single thread.\n")
                num_threads = 1


    # Figure out how many blocks we're going to download
    num_blocks = len(to_download) // block_size
    if len(to_download) % block_size > 0:
        num_blocks += 1
    print(f"Downloading {num_blocks} blocks of ~{block_size} sequences... ")

    # queue will hold results from each download batch.
    queue = mp.Manager().Queue()
    with mp.Pool(num_threads) as pool:

        # This is a bit obscure. Build a list of args to pass to the pool. Each
        # tuple of args matches the args in _entrez_download_thread.
        # all_args has all len(df) reverse blast runs we want to do.
        all_args = []
        for i in range(0,len(to_download),block_size):
            ids = ",".join(to_download[i:(i+block_size)])
            all_args.append((i,ids,num_tries_allowed,wait_time_between_requests,queue))

        # Black magic. pool.imap() runs a function on elements in iterable,
        # filling threads as each job finishes. (Calls _entrez_download_thread
        # on every args tuple in all_args). tqdm gives us a status bar.
        # By wrapping pool.imap iterator in tqdm, we get a status bar that
        # updates as each thread finishes.
        list(tqdm(pool.imap(_entrez_download_thread,all_args),
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
