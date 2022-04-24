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

import re, os, warnings
import multiprocessing as mp

import sys

def _check_for_match(hit,patterns):
    """
    Check to see if a blast hit matches a pattern in patterns. If no matched
    pattern is called, return None.

    hit: biopython blast hit instance
    patterns: list of compiled regular expressions.

    returns: reverse blast hit call, hit definition, called hit e value, e value for best non-call
    """

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

    return hit_call, hit_def, hit_e_value, next_hit_e_value



def reverse_blast(df,
                  call_dict,
                  ncbi_rev_blast_db=None,
                  local_rev_blast_db=None,
                  ncbi_taxid=None,
                  hitlist_size=50,
                  e_value_cutoff=0.01,
                  gapcosts=(11,1),
                  local_num_threads=-1,
                  **kwargs):
    """
    df: topiary dataframe
    call_dict: dictionary with regular expressions as keys and paralog calls as
               values. Both keys and values must be strings.

               example:
               {"lymphogen antigen 96":"LY96",
                "MD-2":"LY96",
                "lymophogen antigen 86":"LY86",
                "MD-1":"LY86"}

                This would mean hits with 'lymophogen antigen 96' and 'MD-2'
                will map to LY96; hits with 'lymphogen antigen 86' and 'MD-1'
                will map to LY86.

    ncbi_rev_blast_db: database on ncbi against which to blast (incompatible
                       with local_rev_blast_db)
    local_rev_blast_db: local database against which to blast (incompatible with
                        ncbi_rev_blast_db)
    ncbi_taxid: limit search to species specified by taxid for an ncbi search
    histlist_size: number of hits to look at
    e_value_cutoff: minimum allowable e value for a hit
    gapcosts: gap costs (must be length 2 tuple of ints)
    local_num_threads: number of threads to use. if -1, use all available.
    kwargs: extra keyword arguments are passed directly to biopython
            NcbiblastXXXCommandline (for local blast) or qblast (for remote
            blast). These take precedence over anything specified above
            (hitlist_size, for example).
    """

    if type(df) is not pd.DataFrame:
        err = "\ndf should be a topiary dataframe\n\n"
        raise ValueError(err)

    # Make sure this is a clean topiary dataframe
    df = util.check_topiary_dataframe(df)

    # Create list of all sequences in dataframe
    sequence_list = []
    for i in range(len(df)):

        # Get sequence
        idx = df.index[i]
        s = df.loc[idx,"sequence"]

        # Try to get start/end (only blasts subset of sequences). If start and
        # end are not defined in this dataframe, take whole sequence
        try:
            a = df.loc[idx,"start"]
            b = df.loc[idx,"end"]
        except KeyError:
            a = 0
            b = None

        # Record this sequence
        sequence_list.append(s[a:b])

    # Check the call_dict argument for sanity
    good_call_dict = True
    if type(call_dict) is not dict:
        good_call_dict = False
    else:
        for k in call_dict:
            if (type(k) is not str) or (call_dict[k] is not str):
                good_call_dict = False
                break
    if not good_call_dict:
        err = "\ncall_dict should be a dictionary keying patterns to look for\n"
        err += "in the blast hits to the paralog call. Both keys and values\n"
        err += "must be strings.\n\n"
        raise ValueError(err)

    # Compile patterns to look for in the blast hits
    patterns = []
    for k in call_dict:
        if type(k) is not str:
            err = "call_dict should be a diti"

        patterns.append((re.compile(k,re.IGNORECASE),call_dict[k]))

    # Figure out what we're blasting against
    if ncbi_rev_blast_db is None and local_rev_blast_db is None:
        err = "\nPlease specificy either ncbi_rev_blast_db OR local_rev_blast_db\n\n"
        raise ValueError(err)

    if ncbi_rev_blast_db is not None and local_rev_blast_db is not None:
        err = "\nPlease specificy either ncbi_rev_blast_db OR\n"
        err += "local_rev_blast_db, but not both.\n\n"
        raise ValueError(err)

    # NCBI blast
    if ncbi_rev_blast_db:
        hit_dfs = ncbi.ncbi_blast(sequence_list,
                                  db=ncbi_rev_blast_db,
                                  taxid=ncbi_taxid,
                                  blast_program="blastp",
                                  hitlist_size=hitlist_size,
                                  e_value_cutoff=e_value_cutoff,
                                  gapcosts=gapcosts,
                                  **kwargs)

    # Local blast
    else:
        hit_dfs = ncbi.local_blast(sequence_list,
                                   db=local_rev_blast_db,
                                   blast_program="blastp",
                                   hitlist_size=hitlist_size,
                                   e_value_cutoff=e_value_cutoff,
                                   gapcosts=gapcosts,
                                   num_threads=local_num_threads,
                                   **kwargs)

    # hit_dfs has a dataframe of blast hits for each sequence in the input
    # topiary dataframe.

    results = []
    for h in hits:
        results.append(_check_for_match(h,patterns))

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
