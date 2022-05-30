__description__ = \
"""
Remove redundancy for datasets in a semi-intelligent way.
"""
__author__ = "Michael J. Harms"
__date__ = "2021-04-08"

from topiary import _arg_processors

import pandas as pd
import numpy as np

from Bio import pairwise2

from tqdm.auto import tqdm

import re

# Columns to check in order
_EXPECTED_COLUMNS = ["structure","low_quality","partial","predicted",
                     "precursor","hypothetical","isoform","diff_from_median",
                     "length"]
_LENGTH_COLUMN = -1

class DummyTqdm():
    def __init__(self):
        pass
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        pass

    def update(self,value):
        pass


def _get_quality_scores(row,key_species={}):
    """
    Get stats in order of importance (see remove_redundancy doc string).

    row: row from dataframe built by topiary
    key_species: dictionary of key species to prefer to others. only uses
                 keys for fast look up and ignores values

    returns float array for the sequence
    """

    # See if this is a key species
    try:
        key_species[row.species]
        values = [0]
    except KeyError:
        values = [1]

    # Add values for expected columns
    expected = list(row[_EXPECTED_COLUMNS])

    # Flip length column
    expected[_LENGTH_COLUMN] = 1/expected[_LENGTH_COLUMN]

    # Combine key species call with expected columns
    values.extend(expected)

    return np.array(values,dtype=float)


def _compare_seqs(A_seq,B_seq,A_qual,B_qual,cutoff,discard_key=False):
    """
    Compare sequence A and B based on alignment. If the sequences are
    similar within cutoff, compare A_stats and B_stats and take the sequence
    with the lower score. Scores have left-right priority.  Will select
    sequence with the first element with a lower score. If sequences are
    similar within cutoff and have equal scores, choose A.

    A_seq: sequence A
    B_seq: sequence B
    A_qual: quality scores for A
    B_qual: quality scores for B
    cutoff: cutoff for sequence comparison (~seq identity. between 0 and 1)
    discard_key: whether or not to discard key species, regardless of their
                 qualities.

    returns bool, bool

    True, True: keep both
    True, False: keep A
    False, True: keep B
    """

    # Get a normalized score: matches/len(shortest)
    score = pairwise2.align.globalxx(A_seq,B_seq,score_only=True)
    norm = score/np.min((len(A_seq),len(B_seq)))

    # If sequence similarity is less than the cutoff, keep both
    if norm <= cutoff:
        return True, True

    # If sequence similarity is greater than the cutoff, select one.
    else:

        # If we are not discarding key sequences and both sequences are
        # from key species, automatically keep both.
        if not discard_key:
            if A_qual[0] == 0 and B_qual[0] == 0:
                return True, True

        # Compare two vectors. Identify first element that differs.

        # Return
        #   True, False if A is smaller
        #   False, True if B is smaller
        #   True, False if equal

        comp = np.zeros(A_qual.shape[0],dtype=np.int8)
        comp[B_qual > A_qual] = 1
        comp[B_qual < A_qual] = -1
        diffs = np.nonzero(comp)[0]

        # No difference, keep A arbitrarily
        if diffs.shape[0] == 0:
            return True, False

        # B > A at first difference, keep A
        elif comp[diffs[0]] > 0:
            return True, False

        # B < A at first difference, keep B
        else:
            return False, True


def remove_redundancy(df,cutoff=0.95,key_species=[],status_bar=True):
    """
    De-duplicate sequences according to cutoff and semi-intelligent heuristic
    criteria.

    Returns a copy of df in which "keep" is set to False for duplicates.

    This intelligently chooses between the two sequences. It favors sequences
    according to specific criteria (in this order of importance):

        1. whether sequence is from a key species
        2. whether it's annotated as structure
        3. how different the sequence length is from the median length
        4. whether it's annotated as low quality
        5. whether it's annotated as partial
        6. whether it's annotated as precursor
        7. whether it's annotated as hypothetical
        8. whether it's annotated as isoform
        9. sequence length (preferring longer)

    df: data frame with sequences.
    cutoff: %identity cutoff for combining removing sequences (0-1)
    key_species: list of key species to prefer.
    """

    # Process arguments
    df = _arg_processors.process_topiary_dataframe(df)
    cutoff = _arg_processors.process_float(cutoff,
                                           "cutoff",
                                           minimum_allowed=0,
                                           maximum_allowed=1)
    key_species = _arg_processors.process_iter(key_species,
                                               "key_species",
                                               required_value_type=str,
                                               is_not_type=[str,dict])
    status_bar = _arg_processors.process_bool(status_bar,"status_bar")

    # Encode key species as a dictionary for fast look up
    key_species = dict([(k,None) for k in key_species])

    starting_keep_number = np.sum(df.keep)

    # If not more than one seq, don't do anything
    if len(df) < 2:
        return df

    # Make sure the dataframe has the columns needed for this comparison. If
    # the dataframe
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

    # Figure out how different each sequence is from the median length.  We
    # want to favor sequences that are closer to the median length than
    # otherwise.
    lengths = df.loc[df.keep,"length"]
    counts, lengths = np.histogram(lengths,bins=int(np.round(2*np.sqrt(len(lengths)),0)))
    median_length = lengths[np.argmax(counts)]
    df["diff_from_median"] = np.abs(df.length - median_length)

    # Get quality scores for each sequence
    quality_scores = []
    for i in range(len(df)):
        quality_scores.append(_get_quality_scores(df.iloc[i,:],key_species))

    unique_species = np.unique(df.species)

    total_calcs = 0
    for s in unique_species:
        a = np.sum(df.species == s)
        total_calcs += a*(a - 1)//2


    if total_calcs > 0:

        print("Removing redundant sequences within species.",flush=True)

        # If a status bar is requested, make one. Otherwise, use dummy version
        if status_bar:
            status = tqdm(total=total_calcs)
        else:
            status = DummyTqdm()

        with status as pbar:

            for s in unique_species:

                # species indexes are iloc row indexes corresponding to species of
                # interest.
                species_mask = df.species == s
                species_rows = np.arange(species_mask.shape[0],dtype=np.uint)[species_mask]
                num_this_species = np.sum(species_mask)
                for x in range(num_this_species):

                    # Get species index (i). If it's already set to keep = False,
                    # don't compare to other sequences.
                    i = df.index[species_rows[x]]
                    if not df.loc[i,"keep"]:
                        continue

                    # Get sequence and quality score for A
                    A_seq = df.loc[i,"sequence"]
                    A_qual = quality_scores[species_rows[x]]

                    # Loop over other sequences in this species
                    for y in range(x+1,num_this_species):

                        # Get species index (j). If it's already set to keep = False,
                        # don't compare to other sequences.
                        j = df.index[species_rows[y]]
                        if not df.loc[j,"keep"]:
                            continue

                        # Get sequence and quality score for B
                        B_seq = df.loc[j,"sequence"]
                        B_qual = quality_scores[species_rows[y]]

                        # Decide which sequence to keep (or both). Discard sequences
                        # even if they are from key species--removing redundancy within
                        # a given species.
                        A_bool, B_bool = _compare_seqs(A_seq,B_seq,A_qual,B_qual,cutoff,
                                                       discard_key=True)

                        # Update keep for each sequence
                        df.loc[i,"keep"] = A_bool
                        df.loc[j,"keep"] = B_bool

                        # If we got rid of A, break out of this loop.  Do not need to
                        # compare to A any more.
                        if not A_bool:
                            break

                pbar.update(num_this_species*(num_this_species - 1)//2)


    N = len(df)
    total_calcs = N*(N-1)//2

    if total_calcs > 0:

        print("Removing redundant sequences, all-on-all.",flush=True)

        # If a status bar is requested, make one. Otherwise, use dummy version
        if status_bar:
            status = tqdm(total=total_calcs)
        else:
            status = DummyTqdm()

        with status as pbar:

            counter = 1
            for x in range(len(df)):

                i = df.index[x]

                # If we've already decided not to keep i, don't even look at it
                if not df.loc[i,"keep"]:
                    pbar.update(N - counter)
                    counter += 1
                    continue

                # Get sequence of sequence and quality scores for i
                A_seq = df.loc[i,"sequence"]
                A_qual = quality_scores[x]

                for y in range(i+1,len(df)):

                    j = df.index[y]

                    # If we've already decided not to keep j, don't even look at it
                    if not df.loc[j,"keep"]:
                        continue

                    # Get sequence of sequence and quality scores for j
                    B_seq = df.loc[j,"sequence"]
                    B_qual = quality_scores[y]

                    # Decide which sequence to keep (or both). Do not discard any
                    # sequence from a key species.
                    A_bool, B_bool = _compare_seqs(A_seq,B_seq,A_qual,B_qual,cutoff,
                                                   discard_key=False)

                    # Update keep for each sequence
                    df.loc[i,"keep"] = A_bool
                    df.loc[j,"keep"] = B_bool

                    # If we got rid of A, break out of this loop.  Do not need to
                    # compare to A any more.
                    if not A_bool:
                        break

                pbar.update(N - counter)
                counter += 1

    final_keep_number = np.sum(df.keep)
    print(f"Reduced {starting_keep_number} --> {final_keep_number} sequences.",flush=True)
    print("Done.",flush=True)

    return df

def find_cutoff(df,
                min_cutoff=0.85,
                max_cutoff=0.99,
                try_n_values=5,
                target_number=500,
                sample_size=100,
                key_species=[]):
    """
    Find an identity cutoff that yields approximately target_number sequences.

    Parameters
    ----------
        min_cutoff: minimum identity cutoff to use
        max_cutoff: maximum identity cutoff to use
        try_n_values: try this many different cutoffs between min and max cutoff
        target_number: find a cutoff that gets approximately this number of
                       sequences
        sample_size: grab this number of sequences from the dataframe to find
                     the cutoff
        key_species: key species to pass to remove_redundancy

    Return
    ------
        cutoff (float) that yields approximately target_number sequences from df
    """

    df = _arg_processors.process_topiary_dataframe(df)
    min_cutoff = _arg_processors.process_float(min_cutoff,
                                               "min_cutoff",
                                               minimum_allowed=0,
                                               maximum_allowed=1)
    max_cutoff = _arg_processors.process_float(max_cutoff,
                                               "max_cutoff",
                                               minimum_allowed=0,
                                               maximum_allowed=1)
    try_n_values = _arg_processors.process_int(try_n_values,
                                               "try_n_values",
                                               minimum_allowed=2)
    target_number = _arg_processors.process_int(target_number,
                                                "target_number",
                                                minimum_allowed=1)
    sample_size = _arg_processors.process_int(sample_size,
                                              "sample_size",
                                              minimum_allowed=1)
    key_species = _arg_processors.process_iter(key_species,
                                               "key_species",
                                               required_value_type=str,
                                               is_not_type=[str,dict])

    # Deal with cutoffs
    if min_cutoff > max_cutoff:
        err = "\nmin_cutoff must be less than or equal to max_cutoff\n\n"
        raise ValueError(err)

    # If min_cutoff and max_cutoff are the same, just return it
    if min_cutoff == max_cutoff:
        return min_cutoff
    else:
        step = (max_cutoff - min_cutoff)/(try_n_values - 1)
        cutoffs = np.array([min_cutoff + step*i for i in range(try_n_values)])

    # Make sure sample size is sane
    if sample_size > len(df.index):
        sample_size = len(df.index)

    # Create reduced dataframe to play with
    to_take = np.random.choice(df.index,sample_size,replace=False)
    small_df = df.loc[to_take,:]

    predicted_keep = []
    with tqdm(total=len(cutoffs)) as pbar:
        for c in cutoffs:
            lower_df = remove_redundancy(small_df,
                                         cutoff=c,
                                         key_species=key_species,
                                         status_bar=False)
            fx_kept = np.sum(lower_df.keep)/np.sum(small_df.keep)
            predicted_keep.append(fx_kept*np.sum(df.keep))

            pbar.update(1)

    predicted_keep = np.array(predicted_keep)

    # If target_number is smaller than what we got with min cutoff, return min
    # cutoff. That's as low as we can go.
    if predicted_keep[0] >= target_number:
        return min_cutoff

    # If our target number is bigger than our predicted number with the max
    # cutoff, return the maxium cutoff
    if predicted_keep[-1] <= target_number:
        return max_cutoff

    # Get the predicted_keep indexe just less than the target
    try:
        bottom = np.where(predicted_keep < target_number)[0][-1]
    except IndexError:
        bottom = None

    # Get the predicted_keep index just more than the target
    try:
        top = np.where(predicted_keep >= target_number)[0][0]
    except IndexError:
        top = None

    # Note: top and bottom should never be None at the same time because we
    # first validated that we were not entirely above the target number above.

    # If we failed to find bottom, just return top
    if bottom is None:
        return cutoffs[top]

    # If we failed to find top, just return bottom
    if top is None:
        return cutoffs[bottom]

    # Interpolate between the values just above and below the top to get cutoff
    rise = (cutoffs[top] - cutoffs[bottom])
    run = (predicted_keep[top] - predicted_keep[bottom])
    slope = rise/run
    intercept = cutoffs[top] - slope*predicted_keep[top]

    return target_number*slope + intercept
