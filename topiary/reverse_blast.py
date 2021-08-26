__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
Reverse blast sequence datasets.
"""

from . import ncbi

import pandas as pd
import numpy as np

from tqdm.auto import tqdm

import re, os
import multiprocessing as mp

import sys

def _reverse_blast_thread(args):
    """
    Run reverse blast on a thread. Should only be called via reverse_blast.

    takes args which are interpreted as (df,i,rev_blast_db,patterns,queue)

    df: expects df has "sequence", "start", and "end" columns. Will return a
        copy of the df with "rev_hit" and "paralog" columns, corresponding
        to top hit title and call based on rev_blast_dict. It will also
        update "keep" to be False for any paralog = None sequences.
    i: iloc index to grab
    rev_blast_db: reverse blast database
    patterns: list of regular expression patterns to match sequences
    queue: multiprocessing queue for storing results
    """

    # parse args
    df = args[0]
    i = args[1]
    rev_blast_db = args[2]
    patterns = args[3]
    queue = args[4]

    index = df.index[i]

    s = df.loc[index,"sequence"]
    a = df.loc[index,"start"]
    b = df.loc[index,"end"]

    seq = s[a:b]
    hit = ncbi.local_blast(seq,db=rev_blast_db,hitlist_size=50)

    hit_call = None
    hit_def = None
    hit_e_value = None
    next_hit_e_value = None
    try:

        # Look at the top hit. See if it matches one of the patterns.
        hit_def = hit.loc[0,"hit_def"]
        for p in patterns:
            if p[0].search(hit_def):
                hit_call = p[1]
                break

        # If we made a hit call...
        if hit_call is not None:

            # Get e-value for the top hit
            hit_e_value = hit.loc[0,"e_value"]

            # Go through the next hits sequentially, looking for the first
            # hit that does not match the search pattern.  Record that
            # match e value. The ratio between hit_e_value and
            # next_hit_e_value tells us whether we are confident in the
            # reverse blast.
            for j in range(1,len(hit)):
                try:
                    next_hit_def = hit.loc[j,"hit_def"]
                    found_match = False
                    for p in patterns:
                        if p[0].search(next_hit_def):
                            found_match = True
                            break

                    # If we did not find a match, record the e-value
                    if not found_match:
                        next_hit_e_value = hit.loc[j,"e_value"]
                        break

                except KeyError:
                    # No more hits!
                    break

    except KeyError:
        pass

    # update queue
    queue.put((i,hit_def,hit_call,hit_e_value,next_hit_e_value))


def reverse_blast(df,call_dict=None,rev_blast_db="GRCh38",num_threads=-1):
    """
    df: expects df has "sequence", "start", and "end" columns. Will return a
        copy of the df with new columns:

        rev_hit: top reverse blast hit
        paralog: paralog, if pattern from call_dict matched top hit
        rev_e_value: e-value for top reverse hit
        next_rev_e_value: e-value for best reverse hit that does *not* match a
                          pattern.

        It will also update "keep" to be False for any paralog = None sequences.
    call_dict: dictionary with regular expressions as keys and calls as values.

               example:
               {"lymphogen antigen 96":"LY96",
                "MD-2":"LY96",
                "lymophogen antigen 86":"LY86",
                "MD-1":"LY86"}

    rev_blast_db: pointer to local blast database for reverse blasting.
    num_threads: number of threads to use. if -1, use all available.
    """

    print("Performing reverse blast...")

    patterns = []
    if call_dict is not None:
        for k in call_dict:
            patterns.append((re.compile(k,re.IGNORECASE),call_dict[k]))

    # Figure out number of threads to use
    if num_threads < 0:
        try:
            num_threads = mp.cpu_count()
        except NotImplementedError:
            num_threads = os.cpu_count()
            if num_threads is None:
                warning.warning("Could not determine number of cpus. Using single thread.\n")
                num_threads = 1

    # queue will hold results from each run.
    queue = mp.Manager().Queue()
    with mp.Pool(num_threads) as pool:

        # This is a bit obscure. Build a list of args to pass to the pool. Each
        # tuple of args matches the args in _reverse_blast_thread.
        # all_args has all len(df) reverse blast runs we want to do.
        all_args = [(df,i,rev_blast_db,patterns,queue) for i in range(len(df))]

        # Black magic. pool.imap() runs a function on elements in iterable,
        # filling threads as each job finishes. (Calls _reverse_blast_thread
        # on every args tuple in all_args). tqdm gives us a status bar.
        # By wrapping pool.imap iterator in tqdm, we get a status bar that
        # updates as each thread finishes.
        list(tqdm(pool.imap(_reverse_blast_thread,all_args),total=len(all_args)))

    # Get results out of the queue.
    results = []
    while not queue.empty():
        results.append(queue.get())

    # Sort results
    results.sort()
    rev_hit = [r[1] for r in results]
    paralog = [r[2] for r in results]
    rev_e_value = [r[3] for r in results]
    next_rev_e_value = [r[4] for r in results]

    new_df = df.copy()
    new_df["rev_hit"] = rev_hit
    new_df["paralog"] = paralog
    new_df["rev_e_value"] = rev_e_value
    new_df["next_rev_e_value"] = next_rev_e_value

    # Remove sequences that do not reverse blast from consideration
    mask = np.array([p is None for p in paralog],dtype=np.bool)
    new_df.loc[mask,"keep"] = False

    print("Done.")

    return new_df
