"""
Polish a near final alignment by removing sequences with long insertions or that
are missing sequence.
"""

import topiary
from topiary._private import check
from topiary.quality.alignment import score_alignment

import numpy as np
import pandas as pd

def _get_cutoff(x,avg_bin_contents=None,pct=0.975):
    """
    Get cutoff corresponding to percentile.

    Parameters
    ----------
    x : numpy.ndarray
        float array to test
    avg_bin_contents : int, optional
        create histograms with widths of len(x)/avg_bin_contents and then define
        percentile based on those histograms. See note.
    pct : float, default=0.975
        get value in x corresponding to percentile.

    Notes
    -----
    If avg_bin_contents is specified, use a coarse-grained histogram approach
    to percentile calculation. This tends to be conservative. If a large number
    of values are in a bin with a high value, the function will end up grabbing
    the cutoff corresponding to the first sparsely populated bin, even if this
    means selecting fewer sequences than would strictly correspond to the top
    pct percentile.
    """

    if avg_bin_contents is None:

        # Get cutoff corresponding to the percentile by sorting
        x = np.array(x.copy())
        x.sort()
        idx = int(np.round(len(x)*pct,0))
        return x[idx]

    else:

        check.check_int(avg_bin_contents,"avg_bin_contents",minimum_allowed=1)

        bins = int(np.round(len(x)/avg_bin_contents,0))
        counts, edges = np.histogram(x,bins=bins)
        mids = (edges[1:] - edges[:-1])/2 + edges[:-1]
        cumsum = np.cumsum(counts)/np.sum(counts)

        return np.min(mids[cumsum >= pct])


def polish_alignment(df,
                     realign=True,
                     sparse_column_cutoff=0.80,
                     align_trim=(0,1),
                     fx_sparse_percentile=0.90,
                     sparse_run_percentile=0.90,
                     fx_missing_percentile=0.90):
    """
    Polish a near-final alignment by removing sequences that have long
    insertions or are missing large chunks of the sequence.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe
    realign : bool, default=True
        align after dropping columns
    sparse_column_cutoff : float, default=0.80
        when checking alignment quality, a column is sparse if it has gaps in
        more than sparse_column_cutoff sequences.
    align_trim : tuple, default=(0,1)
        when checking alignment quality, do not score the first and last parts
        of the alignment. Interpreted like a slice, but with percentages.
        (0.0,1.0) would not trim; (0.05,0,98) would trim the first 0.05 off the
        front and the last 0.02 off the back. The default to this function does
        not trim at all.
    fx_sparse_percentile : float, default=0.90
        flag any sequence that is has a fraction sparse above this percentile
        cutoff.
    sparse_run_percentile : float, default=0.90
        flag any sequence that is has total sparse run length above this
        percentile cutoff.
    fx_missing_percentile : float, default=0.90
        flag any sequence that is has a fraction missing above this percentile
        cutoff.

    Notes
    -----
    The alignment is scored using topiary.quality.score_alignment (see that
    docstring for details). Briefly: columns are characterized as either dense
    (most sequences have a non-gap) or sparse (defined as not dense). This call
    is made using the sparse_column_cutoff argument. This function then
    identifies sequences that have many non-gap characters in sparse columns
    overall, sequences with runs of non-gap characters in long runs of sparse
    columns, and sequences that are missing large portions of the dense columns.
    It drops sequences that have BOTH large fx_sparse AND large sparse_run. It
    also drops sequences that have BOTH large fx_missing AND are flagged as
    partial in the original NCBI entry.
    """

    df = check.check_topiary_dataframe(df)
    realign = check.check_int(realign,"realign")
    # sparse_column_cutoff and align_trim checked immediately by score_alignment
    fx_sparse_percentile = check.check_float(fx_sparse_percentile,
                                             "fx_sparse_percentile",
                                             minimum_allowed=0,
                                             maximum_allowed=1)
    sparse_run_percentile = check.check_float(sparse_run_percentile,
                                              "sparse_run_percentile",
                                              minimum_allowed=0,
                                              maximum_allowed=1)
    fx_missing_percentile = check.check_float(fx_missing_percentile,
                                              "fx_missing_percentile",
                                              minimum_allowed=0,
                                              maximum_allowed=1)


    starting_keep = np.sum(df.keep)

    # Score alignment, generating draft alignment if none in the dataframe
    full_df = score_alignment(df,
                              sparse_column_cutoff=sparse_column_cutoff,
                              align_trim=align_trim,
                              silent=False)

    # Look only at kept sequences
    df = full_df.loc[full_df.keep,:]

    # Get worst fx_sparse and sparse_run sequences
    top_fx_sparse = _get_cutoff(df.fx_in_sparse,pct=fx_sparse_percentile)
    top_sparse_run = _get_cutoff(df.sparse_run_length,pct=sparse_run_percentile)

    # If there are no runs, we'll get a top_sparse_run of 0. Set to -1 -- we 
    # should not select on this
    if top_sparse_run == 0:
        top_sparse_run = np.inf

    to_drop_1 = np.logical_and(df.fx_in_sparse >= top_fx_sparse,
                               df.sparse_run_length >= top_sparse_run)

    # Get worst fx_missing and labeled partial
    top_fx_missing = _get_cutoff(df.fx_missing_dense,pct=fx_missing_percentile)
    df.loc[pd.isnull(df.partial),"partial"] = False

    to_drop_2 = np.logical_and(df.fx_missing_dense >= top_fx_missing,
                               df.partial == True)

    # Final mask for dropping
    final_drop = np.logical_or(to_drop_1,to_drop_2)
    new_keep = np.array(df.keep)
    new_keep[final_drop] = False

    # Update full dataframe, making sure we don't drop anything flagged as
    # always_keep
    full_df.loc[full_df.keep,"keep"] = new_keep
    full_df.loc[full_df["always_keep"],"keep"] = True

    current_keep = np.sum(full_df.keep)
    print(f"Reduced {starting_keep} sequences to {current_keep}\n")

    # Realign, if requested
    if realign:
        full_df = topiary.muscle.align(full_df)

    return full_df
