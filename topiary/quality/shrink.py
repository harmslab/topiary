"""
Shrink the sequence database in a rational way.
"""
import topiary
from topiary._private import check
from topiary.quality.taxonomic import get_merge_blocks
from topiary.quality.redundancy import remove_redundancy
from topiary.quality.alignment import score_alignment

import numpy as np
from tqdm.auto import tqdm

def shrink_in_species(df,redundancy_cutoff=0.98):
    """
    Lower sequence redundancy within individual species, ignoring paralog
    annotation.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe
    redundancy_cutoff : float, default=0.98
        merge sequences that have identities greater than or equal to
        redundancy_cutoff.

    Returns
    -------
    df : pandas.DataFrame
        dataframe with keep set to False for redundant sequences
    """

    df = check.check_topiary_dataframe(df)

    redundancy_cutoff = check.check_float(redundancy_cutoff,
                                          "redundancy_cutoff",
                                          minimum_allowed=0,
                                          maximum_allowed=1)

    # Get uid blocks for kept sequences from unique species in dataframe
    kept_df = df.loc[df.keep,:]
    species = np.unique(kept_df.loc[:,"species"])
    merge_blocks = []
    for s in species:
        merge_blocks.append(kept_df.loc[kept_df["species"] == s,"uid"])

    # Figure out how many merge steps we're doing
    num_to_merge = sum([len(m) for m in merge_blocks])
    with tqdm(total=num_to_merge) as pbar:

        uid_to_keep = []

        # Remove redundant sequences within each merge block
        for uid in merge_blocks:

            this_mask = df.loc[:,"uid"].isin(uid)
            this_df = df.loc[this_mask,:]

            this_df = remove_redundancy(this_df,
                                        cutoff=redundancy_cutoff,
                                        silent=True,
                                        discard_key=True)

            uid_to_keep.extend(this_df.loc[this_df.keep,"uid"])
            pbar.update(n=len(uid))

    # Update dataframe keep with merge results
    keep_mask = df.loc[:,"uid"].isin(uid_to_keep)
    df.loc[:,"keep"] = False
    df.loc[keep_mask,"keep"] = True
    if "always_keep" in df:
        df.loc[df.always_keep,"keep"] = True

    return df

def shrink_redundant(df,
                     paralog_column="recip_paralog",
                     weighted_paralog_split=False,
                     merge_block_size=50,
                     redundancy_cutoff=0.98):
    """
    Lower sequence redundancy by sequence identity within taxonomically informed
    blocks.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe
    paralog_column : str, default="recip_paralog"
        column holding preliminary paralog calls.
    weighted_paralog_split : bool, default=False
        when deciding how much of the total budget to assign to each paralog,
        weight the budget by the number of times each paralog is seen. If False,
        (default), split the budget as evenly as possible between the paralogs
        in the dataframe.
    merge_block_size : int, default=50
        create blocks of paralogs merge_block_size out of the species tree
        to do merging based on sequence identity.
    redundancy_cutoff : float, default=0.98
        merge sequences that have identities greater than or equal to
        redundancy_cutoff.

    Returns
    -------
    df : pandas.DataFrame
        dataframe with keep set to False for redundant sequences
    """

    # --------------------------------------------------------------------------
    # Check input arguments

    df = check.check_topiary_dataframe(df)

    if paralog_column not in df.columns:
        err = f"\nparalog_column '{paralog_column}' not in dataframe.\n"
        raise ValueError(err)

    weighted_paralog_split = check.check_bool(weighted_paralog_split,
                                              "weighted_paralog_split")

    merge_block_size = check.check_int(merge_block_size,
                                       "merge_block_size",
                                       minimum_allowed=1)

    redundancy_cutoff = check.check_float(redundancy_cutoff,
                                          "redundancy_cutoff",
                                          minimum_allowed=0,
                                          maximum_allowed=1)

    # --------------------------------------------------------------------------
    # Do merge calculation

    # Get blocks of sequence to merge. This is a dictionary of uid to attempt to
    # merge keyed to paralog
    merge_blocks = get_merge_blocks(df,
                                    target_seq_number=len(df),
                                    paralog_column=paralog_column,
                                    weighted_paralog_split=False,
                                    target_merge_block_size=merge_block_size)

    num_to_merge = sum([len(merge_blocks[p]) for p in merge_blocks])

    with tqdm(total=num_to_merge) as pbar:

        uid_to_keep = []
        for p in merge_blocks:

            # Remove redundant sequences within each merge block
            for m in merge_blocks[p]:

                budget = m[0]
                uid = m[1]

                this_mask = df.loc[:,"uid"].isin(uid)
                this_df = df.loc[this_mask,:]

                this_df = remove_redundancy(this_df,
                                            cutoff=redundancy_cutoff,
                                            silent=True)

                uid_to_keep.extend(this_df.loc[this_df.keep,"uid"])
                pbar.update(n=1)

    # Update dataframe keep with merge results
    keep_mask = df.loc[:,"uid"].isin(uid_to_keep)
    df.loc[:,"keep"] = False
    df.loc[keep_mask,"keep"] = True
    if "always_keep" in df:
        df.loc[df.always_keep,"keep"] = True

    return df


def shrink_aligners(df,
                    target_seq_number,
                    paralog_column="recip_paralog",
                    weighted_paralog_split=False,
                    sparse_column_cutoff=0.80,
                    align_trim=(0.05,0.95)):
    """
    Select sequences that align best within taxonomically informed blocks.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe
    target_seq_number : int
        number of sequences that should be returned in final dataset
    paralog_column : str, default="recip_paralog"
        column holding preliminary paralog calls.
    weighted_paralog_split : bool, default=False
        when deciding how much of the total budget to assign to each paralog,
        weight the budget by the number of times each paralog is seen. If False,
        (default), split the budget as evenly as possible between the paralogs
        in the dataframe.
    sparse_column_cutoff : float, default=0.80
        when checking alignment quality, a column is sparse if it has gaps in
        more than sparse_column_cutoff sequences.
    align_trim : tuple, default=(0.05,0.95)
        when checking alignment quality, do not score the first and last parts
        of the alignment. Interpreted like a slice, but with percentages.
        (0.0,1.0) would not trim; (0.05,0,98) would trim the first 0.05 off the
        front and the last 0.02 off the back.

    Returns
    -------
    df : pandas.DataFrame
        dataframe with keep set to False for redundant sequences
    """

    # --------------------------------------------------------------------------
    # Get taxonomically informed blocks of sequences to merge

    # Get blocks of sequences to merge
    merge_blocks = get_merge_blocks(df,
                                    target_seq_number=target_seq_number,
                                    paralog_column=paralog_column,
                                    weighted_paralog_split=weighted_paralog_split)

    # --------------------------------------------------------------------------
    # Merge within each merge block according to sequences that align the best

    merge_timings = []
    for p in merge_blocks:
        for m in merge_blocks[p]:
            N = len(m[1])
            if N == 0:
                N = 1
            else:
                N = int(np.round(N*np.log(N),0)) + 1
            merge_timings.append(N)


    with tqdm(total=sum(merge_timings)) as pbar:

        uid_to_keep = []
        counter = 0

        # Go through each paralog
        for p in merge_blocks:

            paralog_df = df.loc[df[paralog_column] == p,:]

            # Merge each merge block based on 1) whether sequences are from key
            # species and then 2) which align best.
            for m in merge_blocks[p]:

                budget = m[0]
                uid = m[1]

                this_df = paralog_df.loc[paralog_df.loc[:,"uid"].isin(uid),:]

                # If any key_species sequences in the block, take them.
                in_species_mask = this_df["key_species"]
                this_uid = list(this_df.loc[in_species_mask,"uid"])
                uid_to_keep.extend(this_uid)

                # Update budget and, if budget exhausted, continue to next block
                budget -= len(this_uid)
                if budget < 1:
                    pbar.update(n=merge_timings[counter])
                    counter += 1
                    continue

                # Create new dataframe that now has all key species entries for
                # this paralog, as well as the block sequences
                new_mask = np.logical_or(paralog_df.loc[:,"uid"].isin(uid),
                                         paralog_df["key_species"])
                this_df = paralog_df.loc[new_mask,:]

                # Score alignment
                this_df = score_alignment(this_df,
                                          align_trim=align_trim,
                                          sparse_column_cutoff=sparse_column_cutoff,
                                          silent=True)

                # Drop the key species from the alignment
                this_df = this_df.loc[np.logical_not(this_df["key_species"]),:]

                # Sort by best aligner
                a = np.array(this_df.fx_missing_dense/np.sum(this_df.fx_missing_dense))
                b = np.array(this_df.sparse_run_length/np.sum(this_df.sparse_run_length))

                # Get the uid for the best aligners
                sort_order = np.argsort(a + b)
                this_uid = np.array(this_df.uid.iloc[sort_order])[:budget]
                uid_to_keep.extend(this_uid)

                pbar.update(n=merge_timings[counter])
                counter += 1

    # Update dataframe keep with merge results
    keep_mask = df.loc[:,"uid"].isin(uid_to_keep)
    df.loc[:,"keep"] = False
    df.loc[keep_mask,"keep"] = True
    if "always_keep" in df:
        df.loc[df.always_keep,"keep"] = True

    return df


def shrink_dataset(df,
                   paralog_column="recip_paralog",
                   seqs_per_column=1,
                   max_seq_number=500,
                   redundancy_cutoff=0.90,
                   merge_block_size=50,
                   weighted_paralog_split=False,
                   sparse_column_cutoff=0.80,
                   align_trim=(0.05,0.95)):
    """
    Select a subset of sequences from a topiary dataframe in a rational way.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe
    paralog_column : str, default="recip_paralog"
        column holding preliminary paralog calls.
    seqs_per_column : float, default=1
        aim to have this number of sequences per column in the key species
        sequences. (For example, if the key sequence is 100 amino acids long,
        seqs_per_column=1 would aim for 100 sequences; 2 would aim for 200
        sequences).
    max_seq_number : int, default=500
        maximum number of sequences to get, regardless of seqs_per_column and
        key sequence length.
    redundancy_cutoff : float, default=0.98
        merge sequences from closely related species with sequence identity
        above cutoff.
    merge_block_size : int, default=50
        create blocks of paralogs merge_block_size out of the species tree
        to do merging based on sequence identity.
    weighted_paralog_split : bool, default=False
        when deciding how much of the total budget to assign to each paralog,
        weight the budget by the number of times each paralog is seen. If False,
        (default), split the budget as evenly as possible between the paralogs
        in the dataframe.
    sparse_column_cutoff : float, default=0.80
        when checking alignment quality, a column is sparse if it has gaps in
        more than sparse_column_cutoff sequences.
    align_trim : tuple, default=(0.05,0.95)
        when checking alignment quality, do not score the first and last parts
        of the alignment. Interpreted like a slice, but with percentages.
        (0.0,1.0) would not trim; (0.05,0,98) would trim the first 0.05 off the
        front and the last 0.02 off the back.
    """

    # --------------------------------------------------------------------------
    # Check input arguments

    df = check.check_topiary_dataframe(df)

    try:
        df.loc[:,paralog_column]
    except KeyError:
        err = f"\nparalog_column '{paralog_column}' not found in dataframe.\n\n"
        raise ValueError(err)

    seqs_per_column = check.check_float(seqs_per_column,
                                        "seqs_per_column",
                                        minimum_allowed=0,
                                        minimum_inclusive=False)

    max_seq_number = check.check_int(max_seq_number,
                                     "max_seq_number",
                                     minimum_allowed=1)

    redundancy_cutoff = check.check_float(redundancy_cutoff,
                                          "redundancy_cutoff",
                                          minimum_allowed=0.0,
                                          maximum_allowed=1.0)

    merge_block_size = check.check_int(merge_block_size,
                                       "merge_block_size",
                                       minimum_allowed=1)


    weighted_paralog_split = check.check_bool(weighted_paralog_split,
                                              "weighted_paralog_split")

    # --------------------------------------------------------------------------
    # Get information for doing redundancy reduction

    # Drop unkept columns
    df = df.loc[df.keep,:]

    # Get key species
    key_species = list(set(df.loc[df.key_species,"species"]))
    key_species.sort()

    # Set target seq number based on sequence length
    key_sequences = df.loc[df.key_species,"sequence"]
    target_seq_number = np.max([len(s) for s in key_sequences])*seqs_per_column
    if target_seq_number > max_seq_number:
        target_seq_number = max_seq_number

    target_seq_number = int(round(target_seq_number*1.1,0))

    print(f"Will build final alignment with ~{target_seq_number} sequences.\n",
          flush=True)

    # How many sequences we start with
    starting_keep = np.sum(df.keep)

    print(f"Number of sequences: {starting_keep}",flush=True)

    # ----------------------------------------------
    # Remove nearly identical sequences with species

    print("Removing redundant sequences within species.",flush=True)

    df = shrink_in_species(df,redundancy_cutoff=redundancy_cutoff)

    # Sanity check -- make sure something is left.
    if np.sum(df.keep) == 0:
        err = "redundancy pass removed all sequences!\n"
        raise ValueError(err)

    current_keep = np.sum(df.keep)
    print(f"Number of sequences: {current_keep}",flush=True)

    # --------------------------------------------------------------------------
    # Remove nearly identical sequences from related species for each paralog

    print("Lowering redundancy based on sequence similarity.",flush=True)

    df = shrink_redundant(df,
                          paralog_column=paralog_column,
                          weighted_paralog_split=False,
                          merge_block_size=merge_block_size,
                          redundancy_cutoff=redundancy_cutoff)

    # Sanity check -- make sure something is left.
    if np.sum(df.keep) == 0:
        err = "redundancy pass removed all sequences!\n"
        raise ValueError(err)

    current_keep = np.sum(df.keep)
    print(f"Number of sequences: {current_keep}",flush=True)

    # --------------------------------------------------------------------------
    # Select best aligners within merge blocks

    print("Selecting best aligning sequences from across species tree.",
          flush=True)

    df = shrink_aligners(df,
                         target_seq_number=target_seq_number,
                         paralog_column=paralog_column,
                         weighted_paralog_split=weighted_paralog_split,
                         sparse_column_cutoff=sparse_column_cutoff,
                         align_trim=align_trim)

    current_keep = np.sum(df.keep)
    print(f"Number of sequences: {current_keep}",flush=True)

    print(f"Reduced {starting_keep} sequences to {current_keep}\n")

    return df
