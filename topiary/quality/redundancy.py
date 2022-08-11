"""
Remove redundancy for datasets in a semi-intelligent way.
"""

import topiary
from topiary._private import check
from topiary._private import threads
from topiary._private import interface

from Bio import pairwise2

import pandas as pd
import numpy as np
from tqdm.auto import tqdm

import os

# Columns to check in order
_EXPECTED_COLUMNS = ["low_quality","partial","predicted",
                     "precursor","hypothetical","isoform","structure"]

def _get_quality_scores(row,
                        key_species={},
                        target_length=None,
                        pct_length_cutoff=None):
    """
    Get stats in order of importance for a specific sequence. (see
    remove_redundancy doc string).

    Parameters
    ----------
    row : pandas.Series
        row from dataframe built by topiary
    key_species : dict
        dictionary of key species to prefer to others. only uses keys for fast
        look up and ignores values
    target_length : int, optional
        target_length for sequence. Must be specified with pct_length_cutoff
    pct_length_cutoff : float, optional
        if length is within pct_length_cutoff of target_length, give good
        quality score. otherwise, bad.

    Returns
    -------
    scores : numpy.ndarray
        float array for the sequence
    """

    try:
        if row.always_keep:
            values = [0]
        else:
            values = [1]
    except AttributeError:
        values = [1]

    # See if this is a key species
    try:
        key_species[row.species]
        values.append(0)
    except KeyError:
        values.append(1)

    # Put in target length information if requested
    if target_length is not None:
        if pct_length_cutoff is None:
            err = "\n\npct_length_cutoff must be specified if target_length is\n"
            err += "specified.\n"
            raise ValueError(err)

        diff = np.abs(len(row.sequence) - target_length)
        pct_diff = diff/target_length
        if pct_diff < pct_length_cutoff:
            values.append(0)
        else:
            values.append(1)

    # Add values for expected columns
    values.extend(list(row[_EXPECTED_COLUMNS]))

    # Flip length column
    values.append(1/len(row.sequence))

    return np.array(values,dtype=float)

def _construct_args(sequence_array,
                    quality_array,
                    keep_array,
                    cutoff,
                    discard_key,
                    num_threads=-1,
                    progress_bar=True):
    """
    Break sequence_array into a rational number of blocks given the number of
    threads and construct a list of keyword arguments to pass to
    _check_block_redundancy, one entry per block.

    Parameters
    ----------
    sequence_array : numpy.ndarray
        array holding sequences to compare.
    quality_array : numpy.ndarray
        array holding vectors of quality scores, one vector for each sequence
    keep_array : numpy.ndarray
        array holding whether or not to keep each sequence. boolean.
        This array is updated and is the primary output of this function.
    cutoff : float
        cutoff between 0 and 1 indicating the fractional similarity between two
        sequences above which they are considered redundant.
    discard_key : bool
        if discard_key is False, a redundant sequence will be tossed even if it
        is from a key species
    num_threads : int, default=-1
        number of threads to use. if -1 use all available
    progress_bar : bool, default=True
        whether or not to show a progress bar

    Returns
    -------
    kwargs_list : list
        list of dictionaries of keyword arguments to pass to
        _check_block_redundancy, one entry per block
    num_threads : int
        number of threads to use for the calculation
    """

    # Get number of threads from machine
    num_threads = threads.get_num_threads(num_threads)

    # Determine number of threads useful for this problem. It's not worth
    # chopping up a super small set of comparisons
    max_useful_threads = len(sequence_array)//50
    if max_useful_threads < 1:
        max_useful_threads = 1

    # Set number of threads
    if num_threads > max_useful_threads:
        num_threads = max_useful_threads

    num_blocks = (num_threads**2 - num_threads)//2 + num_threads
    block_size = len(sequence_array)//num_threads

    # Calculate window sizes to cover whole L x L array, where L is number
    # of sequences
    windows = np.zeros(num_threads,dtype=int)
    windows[:] = len(sequence_array)//num_threads

    # This call spreads remainder of L/num_threads evenly across first
    # remainder windows
    windows[:(len(sequence_array) % num_threads)] += 1

    # Blocks will allow us to tile over whole redundancy matrix
    kwargs_list = []
    for i in range(num_threads):
        for j in range(i,num_threads):

            i_block = (np.sum(windows[:i]),np.sum(windows[:i+1]))
            j_block = (np.sum(windows[:j]),np.sum(windows[:j+1]))

            kwargs_list.append({"i_block":i_block,
                                "j_block":j_block,
                                "sequence_array":sequence_array,
                                "quality_array":quality_array,
                                "keep_array":keep_array,
                                "cutoff":cutoff,
                                "discard_key":discard_key})

    return kwargs_list, num_threads

def _compare_seqs(A_seq,B_seq,A_qual,B_qual,cutoff,discard_key=False):
    """
    Compare sequence A and B based on alignment. If the sequences are
    similar within cutoff, compare A_stats and B_stats and take the sequence
    with the lower score. Scores have left-right priority.  Will select
    sequence with the first element with a lower score. If sequences are
    similar within cutoff and have equal scores, choose A.

    Parameters
    ----------
    A_seq : array
        sequence A
    B_seq : array
        sequence B
    A_qual : array
        quality scores for A
    B_qual : array
        quality scores for B
    cutoff : float
        cutoff for sequence comparison (~seq identity. between 0 and 1)
    discard_key : bool
        whether or not to discard key species, regardless of their qualities.

    Returns
    -------
    bool, bool
        True, True: keep both
        True, False: keep A
        False, True: keep B
    """

    # Get a normalized score: matches/len(shortest)
    score = pairwise2.align.localxx(A_seq,B_seq,score_only=True)
    norm = score/np.min((len(A_seq),len(B_seq)))

    # If sequence similarity is less than the cutoff, keep both
    if norm <= cutoff:
        return True, True

    # If sequence similarity is greater than the cutoff, select one.
    else:

        # If both always keep, return that we keep both
        if A_qual[0] == 0 and B_qual[0] == 0:
            return True, True

        # If we are not discarding key sequences and both sequences are
        # from key species, automatically keep both.
        if not discard_key:
            if A_qual[1] == 0 and B_qual[1] == 0:
                return True, True

        # Compare two vectors. Identify first element that differs.
        comp = np.zeros(A_qual.shape[0],dtype=np.int8)
        comp[B_qual > A_qual] = 1
        comp[B_qual < A_qual] = -1
        diffs = np.nonzero(comp)[0]

        # No difference, keep A arbitrarily
        if diffs.shape[0] == 0:
            return True, False

        # B > A at first difference, keep A
        elif comp[diffs[0]] == 1:
            return True, False

        # B < A at first difference, keep B
        else:
            return False, True

def _redundancy_thread_function(i_block,
                                j_block,
                                sequence_array,
                                quality_array,
                                keep_array,
                                cutoff,
                                discard_key,
                                lock):
    """
    Check for redundancy within a block of the sequence by sequence matrix.
    Updates keep_array in a thread-safe manner. Generally should be called by
    threads.thread_manager.

    Parameters
    ----------
    i_block : tuple
        tuple holding indexes over which to iterate in i
    j_block : tuple
        tuple holding indexes over which to interate in j
    sequence_array : numpy.ndarray
        array holding sequences to compare.
    quality_array : numpy.ndarray
        array holding vectors of quality scores, one vector for each sequence
    keep_array : numpy.ndarray
        array holding whether or not to keep each sequence. boolean. This array
        is updated and is the primary output of this function.
    cutoff : float
        cutoff between 0 and 1 indicating the fractional similarity between two
        sequences above which they are considered redundant.
    discard_key : bool
        if discard_key is False, a redundant sequence will be tossed even if it
        is from a key species
    lock : multiprocessing.Lock
        lock controls access to keep_array
    """

    # Loop over block in i
    for i in range(i_block[0],i_block[1]):

        # Skip if we already know we're not keeping i
        if not keep_array[i]:
            continue

        # Loop over block in j
        for j in range(j_block[0],j_block[1]):

            # Skip if we're below or at diagonal
            if j <= i:
                continue

            # Skip if we already know we're not keeping j
            if not keep_array[j]:
                continue

            # Make comparison
            i_keep, j_keep = _compare_seqs(sequence_array[i],
                                           sequence_array[j],
                                           quality_array[i],
                                           quality_array[j],
                                           cutoff,
                                           discard_key=discard_key)

            # Both false
            if not i_keep and not j_keep:
                lock.acquire()
                try:
                    keep_array[i] = False
                    keep_array[j] = False
                finally:
                    lock.release()

            # i false
            elif not i_keep:
                lock.acquire()
                try:
                    keep_array[i] = False
                finally:
                    lock.release()

            # j false
            elif not j_keep:
                lock.acquire()
                try:
                    keep_array[j] = False
                finally:
                    lock.release()

            else:
                pass

            # Not keeping i, we don't need to keep going
            if not keep_array[i]:
                break


def remove_redundancy(df,
                      cutoff=0.95,
                      target_length_cutoff=0.25,
                      discard_key=False,
                      silent=False,
                      num_threads=-1):
    """
    Remove redundant sequences according to cutoff and semi-intelligent
    heuristics.

    If two sequences are identical within a sequence cutof, it selects the
    sequence to keep according to the following criteria, in order:

    1. whether sequence is from a key species
    2. whether it is significantly different in length from the median length of
       sequences from key species
    3. whether it's annotated as low quality
    4. whether it's annotated as partial
    5. whether it's annotated as precursor
    6. whether it's annotated as hypothetical
    7. whether it's annotated as isoform
    8. whether it's annotated as structure
    9. sequence length (preferring longer)

    Parameters
    ----------
    df : pandas.DataFrame
        topiary data frame with sequences
    cutoff : float, default=0.95
        %identity cutoff for combining removing sequences (between 0 and 1)
    target_length_cutoff : float, default=0.25
        Give a higher quality to any sequence whose length is within
        target_length_cutoff pct of the median key_species sequence length. To
        disable, set to None.
    discard_key : bool, default=False
        whether or not to discard sequences from key species
    silent : bool, default=False
        whether to print output and use status bars
    only_in_species : bool, default=False
        only reduce redundancy within species; do not compare sequences between
        species
    num_threads : int, default=-1
        number of threads to use. If -1, use all available

    Returns
    -------
    topiary_dataframe : pandas.dataframe
        Copy of df in which "keep" is set to False for redundant sequences
    """

    # Process arguments
    df = check.check_topiary_dataframe(df)
    cutoff = check.check_float(cutoff,
                               "cutoff",
                               minimum_allowed=0,
                               maximum_allowed=1)

    if target_length_cutoff is not None:
        target_length_cutoff = check.check_float(target_length_cutoff,
                                                 "target_length_cutoff",
                                                 minimum_allowed=0.)

    silent = check.check_bool(silent,"silent")

    # If not more than one seq, don't do anything
    if len(df) < 2:
        return df

    # Extract key species from dataframe and encode as dictionary for fast
    # lookup
    try:
        key_mask = df.loc[:,"key_species"] == True
        key_species = np.unique(df.loc[key_mask,"species"])
        key_species = dict([(k,None) for k in key_species])
    except KeyError:
        key_species = {}

    # Make sure the dataframe has the columns needed for this comparison. If
    # the dataframe does not have the column, simply set to False. If the
    # dataframe has the column with an unexpected type, die.
    for e in _EXPECTED_COLUMNS:

        try:
            v = df[e]
            v*1.0

        # Column doesn't exist -- record as False
        except KeyError:
            df[e] = False

        # Column exists but can't be interpreted as float. Throw error.
        except TypeError:
            err = "\nThe remove_redundancy function expects a dataframe with \n"
            err += f"column '{e}' that can be interpreted as a number. This\n"
            err += f"column exists but has datatype '{df.loc[:,e].dtype}'.\n"
            err += "To fix, please rename this column and re-run this function.\n\n"
            raise ValueError(err)

    # Determine median sequence length for the sequences from key species
    key_species_list = list(key_species.keys())
    if len(key_species_list) > 0:

        key_mask = df["species"].isin(key_species_list)
        lengths = np.array([len(s) for s in df.loc[key_mask,"sequence"]],dtype=int)
        hist_counts, hist_lengths = np.histogram(lengths,
                                                 bins=int(np.round(2*np.sqrt(len(lengths)),0)))
        median_length = hist_lengths[np.argmax(hist_counts)]

        target_length = median_length
        pct_length_cutoff = target_length_cutoff

    else:
        target_length = None
        pct_length_cutoff = None

    # Get quality scores for each sequence
    all_quality_array = []
    for i in range(len(df)):
        q = _get_quality_scores(df.iloc[i,:],
                                key_species=key_species,
                                target_length=target_length,
                                pct_length_cutoff=pct_length_cutoff)
        all_quality_array.append(q)

    all_quality_array = np.array(all_quality_array)

    progress_bar = not silent

    # Create sequence, quality, and keep arrays that only include sequences we
    # are already keeping
    sequence_array = np.array(df.loc[df.keep,"sequence"])
    quality_array = all_quality_array[df.keep]
    keep_array = np.ones(len(sequence_array),dtype=bool)

    kwargs_list, num_threads = _construct_args(sequence_array=sequence_array,
                                               quality_array=quality_array,
                                               keep_array=keep_array,
                                               cutoff=cutoff,
                                               discard_key=discard_key,
                                               num_threads=num_threads)

    threads.thread_manager(kwargs_list,
                           _redundancy_thread_function,
                           num_threads,
                           progress_bar=progress_bar,
                           pass_lock=True)

    # Update keep array
    df.loc[df.keep,"keep"] = keep_array

    return df
