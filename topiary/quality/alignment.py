"""
Functions for cleaning up an alignment, removing poorly aligned sequences.
"""

import topiary
from topiary._private import check

import pandas as pd
import numpy as np

import re

# structures to convert amino acid sequences to integers and back
AA = "ACDEFGHIKLMNPQRSTVWY-"
AA_TO_INT = dict([(a,i) for i, a in enumerate(AA)])
INT_TO_AA = list(AA)

def _get_sparse_columns(seqs,sparse_column_cutoff=0.5):
    """
    Get True/False array for whether each column is less than cutoff % gaps.

    Parameters
    ----------
    seqs : numpy.ndarray
        2D numpy array holding sequences represented as integers
    sparse_column_cutoff : float, default=0.5
        column is sparse if it has > sparse_column_cutoff fraction gap

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
    Get run-length encoding of an array.

    Parameters
    ----------
    input_array : numpy.ndarray
        numpy input array

    Returns
    -------
    lengths : numpy.ndarray
        length of each run
    start_positions : numpy.ndarray
        indexes in array where runs start
    run_values : numpy.ndarray
        values of each run

    Notes
    -----
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
    seqs : numpy.ndarray
        2D numpy array holding sequences represented as integers

    Returns
    -------
    seqs : numpy.ndarray
        Copy of seqs array with gaps-only columns removed.
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

    Parameters
    ----------
    seqs : numpy.ndarray
        array of sequences
    cumulative_keep : numpy.ndarray
        array of indexes to be passed on for further analyses
    force_keep : numpy.ndarray
        array of bool forcing specific columns to be kept
    sparse_column_cutoff : float
        a column is called sparse if it has > sparse_column_cutoff gaps
    maximum_sparse_allowed : float
        only keep sequences where less than maximum_sparse_allowed of the
        sequence (fraction) is a sparse column. Between 0 and 1.

    Returns
    -------
    kept : numpy.ndarray
        copy of sequences containing only kept sequences
    cumulative_keep : numpy.ndarray
        copy of cumulative_keep containing only kept sequences
    force_keep : numpy.ndarray
        copy of force_keep containing only kept sequences
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

    Parameters
    ----------
    seqs : numpy.ndarray
        array of sequences
    cumulative_keep : numpy.ndarray
        array of indexes to be passed on for further analyses
    force_keep : numpy.ndarray
        array of bool forcing specific columns to be kept
    sparse_column_cutoff : float
        a column is called sparse if it has > sparse_column_cutoff gaps
    minimum_dense_required: float
        only keep sequences that have sequence covering >= than
        minimum_dense_required of the dense columns. Between 0 and 1.

    Returns
    -------
    kept : numpy.ndarray
        copy of sequences containing only kept sequences
    cumulative_keep : numpy.ndarray
        copy of cumulative_keep containing only kept sequences
    force_keep : numpy.ndarray
        copy of force_keep containing only kept sequences
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

    Parameters
    ----------
    seqs : numpy.ndarray
        array of sequences
    cumulative_keep : numpy.ndarray
        array of indexes to be passed on for further analyses
    force_keep : numpy.ndarray
        array of bool forcing specific columns to be kept
    sparse_column_cutoff : float
        a column is called sparse if it has > sparse_column_cutoff gaps
    long_insertion_length: int
        look at sequences that participate in an insertion with >= than
        long_insertion_length contiguous sparse columns.
    long_insertion_fx_cutoff: float
        remove sequences that have sequence covering >= than
        long_insertion_fx_cutoff of the insertion. Between 0 and 1.

    Returns
    -------
    kept : numpy.ndarray
        copy of sequences containing only kept sequences
    cumulative_keep : numpy.ndarray
        copy of cumulative_keep containing only kept sequences
    force_keep : numpy.ndarray
        copy of force_keep containing only kept sequences
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

def score_alignment(df,
                    alignment_column,
                    sparse_column_cutoff=0.95,
                    align_trim=(0.1,0.9)):
    """
    Calculate alignment quality scores for each sequence in the dataframe. The
    resulting scores are loaded as columns into the dataframe. In all cases, a
    higher score is a worse alignment.

    The three calculated scores are:

    + :code:`fx_in_sparse`: fraction of columns from the sequence that have no
      gap where most other sequences have a gap.
    + :code:`fx_missing_dense`: fraction of columns from the sequence that have
      a gap where most other sequences have no gap.
    + :code:`sparse_run_length`: length of longest insertion in the sequence.
      The insertion does not have to be continuous in sequence. Instead, this
      function finds runs of sparse_columns and then counts how many columns
      within each run come from this sequence. If, for example, there are ten
      sparse columns in a row and a sequence has eight non-gap characters
      over that run of sparse columns, the run length would be eight.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe
    alignment_column : str, default="alignment"
        column in dataframe to look for the alignment
    sparse_column_cutoff : float, default=0.95
        when checking alignment quality, a column is sparse if it has gaps in
        more than sparse_column_cutoff sequences.
    align_trim : tuple, default=(0.05,0.95)
        when checking alignment quality, do not score the first and last parts
        of the alignment. Interpreted like a slice, but with percentages.
        (0.0,1.0) would not trim; (0.05,0,98) would trim the first 0.05 off the
        front and the last 0.02 off the back.

    Returns
    -------
    topiary_dataframe : pandas.DataFrame
        copy of df with three new columns containing scores.
    """

    # check dataframe
    df = check.check_topiary_dataframe(df)

    # Check alignment column validity
    try:
        aln = df.loc[df.keep,alignment_column]
    except KeyError:
        err = f"\ndataframe does not have alignment_column '{alignment_column}'\n\n"
        raise ValueError(err)

    # Get length of sequence column
    try:
        col_lengths = list(set([len(s) for s in aln]))
        if len(col_lengths) != 1:
            err = f"\nall sequences in '{alignment_column}' do not have the same legnth\n\n"
            raise ValueError(err)
        align_length = col_lengths[0]

    except TypeError:
        err = f"\n'{alignment_column}' could not be interpreted as sequences\n\n"
        raise ValueError(err)

    # Check float inputs
    sparse_column_cutoff = check.check_float(sparse_column_cutoff,
                                                         "sparse_column_cutoff",
                                                         minimum_allowed=0,
                                                         maximum_allowed=1)

    align_trim = check.check_iter(align_trim,"align_trim",
                                              minimum_allowed=2,
                                              maximum_allowed=2)

    front_trim = check.check_float(align_trim[0],
                                               "align_trim[0]",
                                               minimum_allowed=0.0,
                                               maximum_allowed=1.0)

    back_trim =  check.check_float(align_trim[1],
                                               "align_trim[1]",
                                               minimum_allowed=0.0,
                                               maximum_allowed=1.0)

    if front_trim >= back_trim:
        err = "\nfront_trim must not overlap with back_trim\n\n"
        raise ValueError(err)

    # Get indexes for front and back trimming
    front_index = int(round(align_length*front_trim,0))
    back_index = int(round(align_length*back_trim,0))
    if front_index == back_index:
        if front_index > 0:
            front_index -= 1
        else:
            back_index += 1

    # Generate an array of sequences as integers.
    seqs = []
    for i in df.index:
        if not df.loc[i,"keep"]:
            continue
        this_seq = re.sub(f"[^{AA}]","-",df.loc[i,alignment_column][front_index:back_index])
        seqs.append([AA_TO_INT[c] for c in list(this_seq)])

    seqs = np.array(seqs,dtype=int)

    # Drop gaps only columns
    seqs = _drop_gaps_only(seqs)

    sparse_columns = _get_sparse_columns(seqs,sparse_column_cutoff)
    non_gap_sparse = np.logical_and(seqs != 20,sparse_columns)
    fx_in_sparse = np.sum(non_gap_sparse,axis=1)/seqs.shape[1]

    dense_columns = np.logical_not(sparse_columns)
    gap_dense = np.logical_and(seqs == 20,dense_columns)
    fx_missing_dense = np.sum(gap_dense,axis=1)/np.sum(dense_columns)

    sparse_run_length = np.zeros(len(seqs))
    run_lengths, start_positions, values = _rle(dense_columns)
    for i in range(len(run_lengths)):

        # If the this is a run of sparse columns
        if values[i] == False:

            # Grab indexes for start and stop of sparse column run
            s = start_positions[i]
            e = s + run_lengths[i]

            # Go through every sequence
            for j in range(len(seqs)):
                in_sparse_run = np.sum(seqs[j,s:e] != 20)
                if in_sparse_run > sparse_run_length[j]:
                    sparse_run_length[j] = in_sparse_run


    # Load quality data into the dataframe
    df["fx_in_sparse"] = np.nan
    df.loc[df.keep,"fx_in_sparse"] = fx_in_sparse
    df["fx_missing_dense"] = np.nan
    df.loc[df.keep,"fx_missing_dense"] = fx_missing_dense
    df["sparse_run_length"] = np.nan
    df.loc[df.keep,"sparse_run_length"] = sparse_run_length

    return df

def clean_alignment(df,
                    alignment_column="alignment",
                    key_species=[],
                    sparse_column_cutoff=0.95,
                    maximum_sparse_allowed=0.025,
                    minimum_dense_required=0.95,
                    long_insertion_length=8,
                    long_insertion_fx_cutoff=0.6):
    """
    Remove poorly aligned sequences from a topiary dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe with alignment.
    alignment_column : str, default="alignment"
        column in dataframe to look for the alignment
    key_species : list, default=[]
        list of species (binomial names) that will not have sequences removed
        regardless of alignment quality.
    sparse_column_cutoff : float, default=0.95
        when checking alignment quality, a column is sparse if it has gaps in
        more than sparse_column_cutoff sequences.
    maximum_sparse_allowed : float, default=0.025
        only keep sequences where less than maximum_sparse_allowed of the
        sequence (fraction) is a sparse column. Between 0 and 1.
    minimum_dense_required : float, default=0.95
        only keep sequences that have sequence covering >=
        than minimum_dense_required of the dense columns. Between 0 and 1.
    long_insertion_length : int, default=8
        look at sequences that participate in an insertion with >= than
        long_insertion_length contiguous sparse columns.
    long_insertion_fx_cutoff : float, default=0.6
        remove sequences that have sequence covering >= than
        long_insertion_fx_cutof of the insertion. Between 0 and 1.

    Return
    ------
    topiary_dataframe : pandas.DataFrame
        dataframe with cleaned up alignment. sequences removed from the
        alignment have keep=False and their aligned sequence set to pandas.NA.
    """

    # check dataframe
    df = check.check_topiary_dataframe(df)

    # Check alignment column validity
    try:
        df.loc[:,alignment_column]
    except KeyError:
        err = f"\ndataframe does not have alignment_column '{alignment_column}'\n\n"
        raise ValueError(err)


    key_species = check.check_iter(key_species)

    # Check numerical arguments

    sparse_column_cutoff = check.check_float(sparse_column_cutoff,
                                                         "sparse_column_cutoff",
                                                         minimum_allowed=0,
                                                         maximum_allowed=1)

    maximum_sparse_allowed = check.check_float(maximum_sparse_allowed,
                                                           "maximum_sparse_allowed",
                                                           minimum_allowed=0,
                                                           maximum_allowed=1)

    long_insertion_length = check.check_int(long_insertion_length,
                                                        "long_insertion_length",
                                                        minimum_allowed=0)

    long_insertion_fx_cutoff = check.check_float(long_insertion_fx_cutoff,
                                                             "long_insertion_fx_cutoff",
                                                             minimum_allowed=0,
                                                             maximum_allowed=1)

    # Convert alignment into array of integers, only considering sequences we
    # were keeping already
    seqs = []
    force_keep = []

    # If the "always_keep" columns exists, keep it
    try:
        always = np.array(df.loc[:,"always_keep"],dtype=bool)
    except KeyError:
        always = np.zeros(len(df),dtype=bool)

    termini = int(round(len(df.loc[df.index[0],alignment_column])*0.05,0))

    for i in df.index:
        if not df.loc[i,"keep"]:
            continue
        this_seq = re.sub(f"[^{AA}]","-",df.loc[i,alignment_column])
        seqs.append([AA_TO_INT[c] for c in list(this_seq)])

        # If the sequence is from a key species, set force_keep to True for that
        # sequence
        if df.loc[i,"species"] in key_species or always[i]:
            force_keep.append(True)
        else:
            force_keep.append(False)



    seqs = np.array(seqs,dtype=int)
    force_keep = np.array(force_keep,dtype=bool)

    starting_keep_number = np.sum(df.keep)

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
    #df.loc[cumulative_keep,alignment_column] = final_seqs
    no_keep = np.array(list(set(df.index).difference(set(cumulative_keep))))
    df.loc[no_keep,alignment_column] = pd.NA

    final_keep_number = np.sum(df.keep)

    print("Removing excessively gapped sequences.")
    print(f"Reduced {starting_keep_number} --> {final_keep_number} sequences.",flush=True)

    return df
