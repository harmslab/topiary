__description__ = \
"""
Tool for cleaning up alignment, removing sequences with long indels.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-05-27"

import topiary
from topiary import _arg_processors

import pandas as pd
import numpy as np

import re

# structures to convert amino acid sequences to integers and back
AA = "ACDEFGHIKLMNPQRSTVWY-"
AA_TO_INT = dict([(a,i) for i, a in enumerate(AA)])
INT_TO_AA = list(AA)

def _get_sparse_columns(seqs,sparse_column_cutoff=0.9):
    """
    Get True/False array for whether each column is less than cutoff % gaps.

    Parameters
    ----------
        seqs: 2D numpy array holding sequences represented as integers
        sparse_column_cutoff: column is sparse if it has > sparse_column_cutoff
                              fraction gap

    Return
    ------
        True/False numpy array mask
    """

    column_scores = []
    for c in range(seqs.shape[1]):
        column_scores.append(np.sum(seqs[:,c] == 20)/seqs.shape[0])
    column_scores = np.array(column_scores)

    sparse_columns = column_scores >= sparse_column_cutoff

    return sparse_columns

def _rle(input_array):
    """
    Get run-length encoding of an array. For input_array, extract runs of
    identical values.

    Parameters
    ----------
        input_array: numpy input array
        get_only_true: if set to True, return only runs of "True" from
                       input_array.

    Return
    ------
        run_lengths, start_positions, values

    This is a cleaned up version of a solution posted here:
    https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encodi
    """

    N = input_array.shape[0]

    # Get differences in input array
    diffs = input_array[1:] != input_array[:-1]

    # Build array to hold indexes of differences. Stick -1 at front and
    # len(input_array) -1 at back so we count from start and end.
    diff_indexes = np.zeros(np.sum(diffs) + 2,dtype=int)
    diff_indexes[0] = -1
    diff_indexes[-1] = N - 1

    # Get the indexes where the sequences differ from each other.
    diff_indexes[1:-1] = np.where(diffs)[0]

    # Get lengths of jumps between indexes
    run_lengths = np.diff(diff_indexes)

    # Start positions array. Append 0 on front -- first always starts at zero.
    # Don't include last index -- not a start.
    start_positions = np.zeros(run_lengths.shape[0],dtype=int)
    start_positions[1:] = np.cumsum(run_lengths[:-1])

    # These are the values of each of the runs
    values = input_array[diff_indexes[1:]]

    return run_lengths, start_positions, values

def _drop_gaps_only(seqs):
    """
    Drop gap-only columns from an array holding a sequence alignment.

    Parameters
    ----------
        seqs: 2D numpy array holding sequences represented as integers

    Return
    ------
        Copy of array with gaps-only columns removed.
    """

    # Create True/False mask for columns with more than just "-"
    not_just_gaps = []
    for i in range(seqs.shape[1]):
        u = np.unique(seqs[:,i])
        if len(u) == 1 and u[0] == 20:
            not_just_gaps.append(False)
        else:
            not_just_gaps.append(True)

    not_just_gaps = np.array(not_just_gaps,dtype=bool)

    # Whack out columns that are only "-"
    seqs = seqs[:,not_just_gaps]

    return seqs

def _find_too_many_sparse(seqs,
                          cumulative_keep,
                          force_keep,
                          sparse_column_cutoff,
                          maximum_sparse_allowed):
    """
    Filter out sequences where over maximum_sparse_allowed of sequence is in
    sparse columns.
    """

    sparse_columns = _get_sparse_columns(seqs,sparse_column_cutoff)
    non_gap_sparse = np.logical_and(seqs != 20,sparse_columns)
    seq_score = np.sum(non_gap_sparse,axis=1)/seqs.shape[1]

    keep_mask = seq_score <= maximum_sparse_allowed
    keep_mask = np.logical_or(keep_mask,force_keep)
    force_keep = force_keep[keep_mask]
    cumulative_keep = cumulative_keep[keep_mask]
    seqs = seqs[keep_mask]
    seqs = _drop_gaps_only(seqs)

    return seqs, cumulative_keep, force_keep

def _find_too_few_dense(seqs,
                        cumulative_keep,
                        force_keep,
                        sparse_column_cutoff,
                        minimum_dense_required):
    """
    Filter out sequences that have a non gap in less than minimum_dense_required
    fraction of the dense columns.
    """

    sparse_columns = _get_sparse_columns(seqs,sparse_column_cutoff)
    dense_columns = np.logical_not(sparse_columns)
    gap_dense = np.logical_and(seqs != 20,dense_columns)
    seq_score = np.sum(gap_dense,axis=1)/np.sum(dense_columns)

    keep_mask = seq_score >= minimum_dense_required
    keep_mask = np.logical_or(keep_mask,force_keep)
    force_keep = force_keep[keep_mask]
    cumulative_keep = cumulative_keep[keep_mask]
    seqs = seqs[keep_mask]
    seqs = _drop_gaps_only(seqs)

    return seqs, cumulative_keep, force_keep

def _find_long_insertions(seqs,
                          cumulative_keep,
                          force_keep,
                          sparse_column_cutoff,
                          long_insertion_length,
                          long_insertion_fx_cutoff):
    """
    Filter out sequences that are part of long stretches of sparse sequences.
    Any sequence with long_insertion_fx_cutoff of the sparse columns not
    gapped will be removed.
    """

    sparse_columns = _get_sparse_columns(seqs,sparse_column_cutoff)
    dense_columns = np.logical_not(sparse_columns)
    num_dense_columns = np.sum(dense_columns)

    run_lengths, start_positions, values = _rle(dense_columns)

    keep_mask = np.ones(seqs.shape[0],dtype=bool)
    for i in range(len(run_lengths)):

        # If the this is a run of sparse columns
        if values[i] == False:

            # If the run length is less than the long_insertion_length stop
            # considering it
            if run_lengths[i] < long_insertion_length:
                continue

            # Grab indexes for start and stop of sparse column run
            s = start_positions[i]
            e = s + run_lengths[i]

            # Go through every sequence
            for j in range(len(seqs)):

                # If we haven't already tossed this sequence...
                if keep_mask[j]:
                    # Calculate the fraction of non-gap characters for this
                    # sequence in the sparse region. If this is greater than
                    # long_insertion_fx_cutoff, toss the sequence
                    fx_in_sparse = np.sum(seqs[j,s:e] != 20)/run_lengths[i]
                    if fx_in_sparse >= long_insertion_fx_cutoff:
                        keep_mask[j] = False

    keep_mask = np.logical_or(keep_mask,force_keep)
    force_keep = force_keep[keep_mask]
    cumulative_keep = cumulative_keep[keep_mask]
    seqs = seqs[keep_mask]
    seqs = _drop_gaps_only(seqs)

    return seqs, cumulative_keep, force_keep

def clean_alignment(df,
                    alignment_column="alignment",
                    key_species=[],
                    sparse_column_cutoff=0.5,
                    maximum_sparse_allowed=0.025,
                    minimum_dense_required=0.90,
                    long_insertion_length=8,
                    long_insertion_fx_cutoff=0.8):
    """

    Parameters
    ----------
        df: topiary dataframe with alignment.
        alignment_column: column in dataframe to look for the alignment
        key_species: list of species whose sequences will be preserved,
                     regardless of their alignment quality.
        sparse_column_cutoff: a column is called sparse if it has
                              > sparse_column_cutoff gaps
        maximum_sparse_allowed: only keep sequences where less than
                                maximum_sparse_allowed of the sequence (fraction)
                                is a sparse column. Between 0 and 1. Default: 0.025.
        minimum_dense_required: only keep sequences that have sequence covering
                                >= than minimum_dense_required of the
                                dense columns. Between 0 and 1. Default: 0.9.
        long_insertion_length: look at sequences that participate in an insertion
                               with >= than long_insertion_length contiguous
                               sparse columns.
        long_insertion_fx_cutoff: remove sequences that have sequence covering
                                  >= than long_insertion_fx_cutof of the
                                  insertion. Between 0 and 1. Default: 0.8.

    Return
    ------
        dataframe with cleaned up alignment. sequences removed from the alignment
        have keep=False and their aligned sequence set to pd.NA.
    """

    # check dataframe
    df = _arg_processors.process_topiary_dataframe(df)

    # Check alignment column validity
    try:
        df.loc[:,alignment_column]
    except KeyError:
        err = f"\ndataframe does not have alignment_column '{alignment_column}'\n\n"
        raise ValueError(err)


    key_species = _arg_processors.process_iter(key_species)

    # Check numerical arguments

    sparse_column_cutoff = _arg_processors.process_float(sparse_column_cutoff,
                                                         "sparse_column_cutoff",
                                                         minimum_allowed=0,
                                                         maximum_allowed=1)

    maximum_sparse_allowed = _arg_processors.process_float(maximum_sparse_allowed,
                                                           "maximum_sparse_allowed",
                                                           minimum_allowed=0,
                                                           maximum_allowed=1)

    long_insertion_length = _arg_processors.process_int(long_insertion_length,
                                                        "long_insertion_length",
                                                        minimum_allowed=0)

    long_insertion_fx_cutoff = _arg_processors.process_float(long_insertion_fx_cutoff,
                                                             "long_insertion_fx_cutoff",
                                                             minimum_allowed=0,
                                                             maximum_allowed=1)

    # Convert alignment into array of integers, only considering sequences we
    # were keeping already
    seqs = []
    force_keep = []
    for i in df.index:
        if not df.loc[i,"keep"]:
            continue
        this_seq = re.sub(f"[^{AA}]","-",df.loc[i,alignment_column])
        seqs.append([AA_TO_INT[c] for c in list(this_seq)])

        # If the sequence is from a key species, set force_keep to True for that
        # sequence
        if df.loc[i,"species"] in key_species:
            force_keep.append(True)
        else:
            force_keep.append(False)

    seqs = np.array(seqs,dtype=int)
    force_keep = np.array(force_keep,dtype=bool)

    # Will hold indexes to keep after filtering, getting rid of those with
    # keep = False from the start
    cumulative_keep = df.index[df.keep]

    # Remove sequences with long runs of non-gap characters in low-quality
    # regions. (Basically, sequences with long, unique insertions)
    seqs, cumulative_keep, force_keep = _find_too_many_sparse(seqs,
                                                              cumulative_keep,
                                                              force_keep,
                                                              sparse_column_cutoff,
                                                              maximum_sparse_allowed)

    # Remove sequences with that do not cover enough of the dense columns
    # (basically, partial sequences)
    seqs, cumulative_keep, force_keep = _find_too_few_dense(seqs,
                                                            cumulative_keep,
                                                            force_keep,
                                                            sparse_column_cutoff,
                                                            minimum_dense_required)


    # Remove sequences that have long stretches of sequence in runs of sparse
    # columns
    seqs, cumulative_keep, force_keep = _find_long_insertions(seqs,
                                                              cumulative_keep,
                                                              force_keep,
                                                              sparse_column_cutoff,
                                                              long_insertion_length,
                                                              long_insertion_fx_cutoff)



    # Construct final set of sequences as strings
    final_seqs = []
    for i in range(seqs.shape[0]):
        final_seqs.append("".join([INT_TO_AA[s] for s in seqs[i,:]]))

    # Construct final dataframe with cleaned up alignment
    df.loc[:,"keep"] = False
    df.loc[cumulative_keep,"keep"] = True
    df.loc[cumulative_keep,alignment_column] = final_seqs
    no_keep = np.array(list(set(df.index).difference(set(cumulative_keep))))
    df.loc[no_keep,alignment_column] = pd.NA

    return df
