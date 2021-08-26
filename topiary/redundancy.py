__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
Remove redundancy for datasets in a semi-intelligent way.
"""

import pandas as pd
import numpy as np

from Bio import pairwise2

from tqdm.auto import tqdm

import re

def _get_quality_scores(row,key_species={}):
    """
    Get stats in order of importance (see remove_redundancy doc string).

    row: row from dataframe built by topiary
    key_species: dictionary of key species to prefer to others. only uses
                 keys for fast look up and ignores values

    returns float array for the sequence
    """

    try:
        key_species[row.species]
        key_species_score = 0.0
    except KeyError:
        key_species_score = 1.0

    return np.array([key_species_score,
                     row.structure,
                     row.diff_from_median,
                     row.low_quality,
                     row.partial,
                     row.precursor,
                     row.hypothetical,
                     row.isoform,
                     1/row.length],dtype=np.float)

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
    norm = score/min((len(A_seq),len(B_seq)))

    # If sequence similarity is less than the cutoff, keep both
    if norm <= cutoff:
        return True, True

    # If sequence similarity is greater than the cutoff, select one.
    else:

        # If we are not discarding key sequences and both sequences are
        # from key species, automatically keep both.
        if not discard_key:
            if A_qual[0] == 1 and B_qual[0] == 1:
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


def remove_redundancy(df,cutoff=0.95,key_species=[]):
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

    key_species = dict([(k,None) for k in key_species])

    # This will hold output
    new_df = df.copy()

    # If not more than one seq, don't do anything
    if len(df) < 2:
        return new_df

    # Figure out how different each sequence is from the median length.  We
    # want to favor sequences that are closer to the median length than
    # otherwise.
    lengths = df.loc[df.keep,"length"]
    counts, lengths = np.histogram(lengths,bins=np.int(np.round(2*np.sqrt(len(lengths)),0)))
    median_length = lengths[np.argmax(counts)]
    new_df["diff_from_median"] = np.abs(new_df.length - median_length)

    # Get quality scores for each sequence
    quality_scores = []
    for i in range(len(new_df)):
        quality_scores.append(_get_quality_scores(new_df.iloc[i,:],key_species))

    print("Removing redundant sequences within species.")
    unique_species = np.unique(new_df.species)

    total_calcs = 0
    for s in unique_species:
        a = np.sum(new_df.species == s)
        total_calcs += a*(a - 1)//2

    with tqdm(total=total_calcs) as pbar:

        for s in unique_species:

            # species indexes are iloc row indexes corresponding to species of
            # interest.
            species_mask = new_df.species == s
            species_rows = np.arange(species_mask.shape[0],dtype=np.uint)[species_mask]
            num_this_species = np.sum(species_mask)
            for x in range(num_this_species):

                # Get species index (i). If it's already set to keep = False,
                # don't compare to other sequences.
                i = df.index[species_rows[x]]
                if not new_df.loc[i,"keep"]:
                    continue

                # Get sequence and quality score for A
                A_seq = new_df.loc[i,"sequence"]
                A_qual = quality_scores[species_rows[x]]

                # Loop over other sequences in this species
                for y in range(x+1,num_this_species):

                    # Get species index (j). If it's already set to keep = False,
                    # don't compare to other sequences.
                    j = df.index[species_rows[y]]
                    if not new_df.loc[j,"keep"]:
                        continue

                    # Get sequence and quality score for B
                    B_seq = new_df.loc[j,"sequence"]
                    B_qual = quality_scores[species_rows[y]]

                    # Decide which sequence to keep (or both). Discard sequences
                    # even if they are from key species--removing redundancy within
                    # a given species.
                    A_bool, B_bool = _compare_seqs(A_seq,B_seq,A_qual,B_qual,cutoff,
                                                   discard_key=True)

                    # Update keep for each sequence
                    new_df.loc[i,"keep"] = A_bool
                    new_df.loc[j,"keep"] = B_bool

                    # If we got rid of A, break out of this loop.  Do not need to
                    # compare to A any more.
                    if not A_bool:
                        break

            pbar.update(num_this_species*(num_this_species - 1)//2)

    print("Removing redundant sequences, all-on-all.")

    N = len(new_df)
    total_calcs = N*(N-1)//2
    with tqdm(total=total_calcs) as pbar:

        counter = 1
        for x in range(len(new_df)):

            i = new_df.index[x]

            # If we've already decided not to keep i, don't even look at it
            if not new_df.loc[i,"keep"]:
                pbar.update(N - counter)
                counter += 1
                continue

            # Get sequence of sequence and quality scores for i
            A_seq = new_df.loc[i,"sequence"]
            A_qual = quality_scores[x]

            for y in range(i+1,len(new_df)):

                j = new_df.index[y]

                # If we've already decided not to keep j, don't even look at it
                if not new_df.loc[j,"keep"]:
                    continue

                # Get sequence of sequence and quality scores for j
                B_seq = new_df.loc[j,"sequence"]
                B_qual = quality_scores[y]

                # Decide which sequence to keep (or both). Do not discard any
                # sequence from a key species.
                A_bool, B_bool = _compare_seqs(A_seq,B_seq,A_qual,B_qual,cutoff,
                                               discard_key=False)

                # Update keep for each sequence
                new_df.loc[i,"keep"] = A_bool
                new_df.loc[j,"keep"] = B_bool

                # If we got rid of A, break out of this loop.  Do not need to
                # compare to A any more.
                if not A_bool:
                    break

            pbar.update(N - counter)
            counter += 1

    print("Done.")

    return new_df
