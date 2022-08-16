"""
Merge dataframes from multiple BLAST queries.
"""

import topiary
from topiary._private import check
from topiary.ncbi.entrez.sequences import get_sequences
from topiary.ncbi import parse_ncbi_line

import numpy as np
import pandas as pd

def _check_merge(merge_this,merge_into,merge_list):
    """
    Merge the row given by index into the row given by merge_index. If the row
    given by merge_index is itself already merging, merge the index row into
    the existing merge row.

    Parameters
    ----------
        index: index of the row we are going to merge
        merge_index: index of the row we are going to merge INTO
        to_merge: list of sequences to merge. if to_merge[merge_index] is
                  already defined, set merge_index -> to_merge[merge_index]
                  and then do to_merge[index] = merge_index

    Return
    ------
        updated to_merge list
    """

    # If merge_list[merge_into] is not defined, set it to merge_into. If it is
    # defined, set merge_into to merge_list[merge_into]
    if merge_list[merge_into] is None:
        merge_list[merge_into] = merge_into
    else:
        merge_into = merge_list[merge_into]

    # If the sequence we're merging has not been merged, just record.
    if merge_list[merge_this] is None:
        merge_list[merge_this] = merge_into
    else:

        # If the current value in merge_list is *not* merge_into, update all
        # entries that have the current value to be merge_into.
        if merge_list[merge_this] != merge_into:
            curr_value = merge_list[merge_this]
            for i in range(len(merge_list)):
                if merge_list[i] in [curr_value,merge_into]:
                    merge_list[i] = merge_into

    return merge_list

def merge_blast_df(blast_df_list):
    """
    Merge dataframes from multiple BLAST queries. Merge happens based on
    accession. If different queries returned an overlapping sequence region
    from the same accession, merge them together. The merge can include an
    arbitrary number of sequences. This merge will:

    1. Take the longest contiguous hit as the row to keep.
    2. Update the "query" column to include all merged queries. For
       example, if the Human and Troll queries gave overlapping hits from
       subject region 1-100 and 2-101, respectively, the merged query will
       be "Human|1-100;Troll|2-101.
    3. Expand the "subject_start" and "subject_end" numbers to cover the
       maximum extent of all hits. In the example from (1), these would
       become 1 and 101 (human start, troll end).
    4. Set the "e_value" column to be the lowest e_value observed by any
       of the merged hits.

    Parameters
    ----------
    blast_df_list : list
        list of dataframes returned by a blast call with multiple sequence
        queries

    Returns
    -------
    merged : pandas.DataFrame
        single dataframe with duplicate hits between queries merged
    """

    # Make sure blast_df_list is iterable
    try:
        blast_df_list = check.check_iter(blast_df_list,
                                                     "blast_df_list")
    except ValueError:
        err = "\nblast_df_list must be a list of pandas dataframes\n\n"
        raise ValueError(err)

    # Validate input
    required_columns = ["accession","subject_start","subject_end","query","e_value"]
    for d in blast_df_list:

        # Make sure elements are all dataframes
        if type(d) is not type(pd.DataFrame()):
            err = "\nblast_df_list must be a list of pandas dataframes\n\n"
            raise ValueError(err)

        # Make sure the dataframe has the required columns
        for r in required_columns:
            if not r in d.columns:
                err = f"\ndataframe in blast_df_list does not have required column '{r}'.\n"
                err += "Dataframes must have columns:\n"
                for c in required_columns:
                    err += f"    {c}\n"
                err += "\n"
                raise ValueError(err)

    # Don't do anything if we send in a single entry -- just return a
    # copy of the dataframe. (We send out a copy because the merged dataframe
    # will generally not be the same object as the inputs).
    if len(blast_df_list) == 0:
        return blast_df_list[0].copy()

    # concatenate dataframes
    df = pd.concat(blast_df_list,ignore_index=True)

    # Mask will hold whether to keep each row or not
    keep_mask = np.ones(len(df.index),dtype=bool)

    # Go over all accessions
    accessions = np.unique(df.loc[:,"accession"])
    for a in accessions:

        # Grab every row with this accession
        mask = df.accession == a

        # If only one row, move on
        num_entries = np.sum(mask)
        if num_entries == 1:
            continue

        # Get copy of dataframe with only this accession
        this_df = df.loc[mask,:]

        # Compare all entries that have the same accession
        to_merge = [None for _ in range(num_entries)]
        for i in range(num_entries):

            i_start = this_df.subject_start.iloc[i]
            i_end = this_df.subject_end.iloc[i]

            for j in range(i+1,num_entries):

                j_start = this_df.subject_start.iloc[j]
                j_end = this_df.subject_end.iloc[j]

                # If they overlap, merge them
                if not (j_end < i_start or i_end < j_start):
                    to_merge = _check_merge(i,j,to_merge)

        # Get unique overlapping seuqneces
        clusters = list(set(to_merge))
        for cluster in clusters:

            # Don't merge; no overlaps detected
            if cluster is None:
                continue

            # Get all rows to merge, including the one to keep
            merge_idx = np.array([cluster == c for c in to_merge],dtype=bool)
            rows_to_merge = this_df.index[merge_idx]
            keep_mask[rows_to_merge] = False

            # Get the one row to keep
            row_to_keep = this_df.index[cluster]
            keep_mask[row_to_keep] = True

            # Create new start, ends, querie, and e-value for row to keep

            # Get starts, ends, queries, and e-value for the rows to merge
            starts = this_df.loc[rows_to_merge,"subject_start"]
            ends = this_df.loc[rows_to_merge,"subject_end"]
            queries = this_df.loc[rows_to_merge,"query"]
            e_value = this_df.loc[rows_to_merge,"e_value"]

            # Broadest start to end from the merged set
            new_start = np.min(starts)
            new_end = np.max(ends)

            # Lowest e-value from merged set
            new_evalue = np.min(e_value)

            # Create a new query, merging the queries and starts/stops from the
            # merged rows.
            new_query = []
            for i in range(len(queries)):
                new_query.append(f"{queries.iloc[i]}|{starts.iloc[i]}-{ends.iloc[i]}")
            new_query = ";".join(new_query)

            # Update the row to keep with the new query, start, and end
            df.loc[row_to_keep,"query"] = new_query
            df.loc[row_to_keep,"subject_start"] = new_start
            df.loc[row_to_keep,"subject_end"] = new_end
            df.loc[row_to_keep,"e_value"] = new_evalue


    return df.loc[keep_mask,:].reset_index().drop(columns=["index"])


def merge_and_annotate(blast_df_list,blast_source_list=None):
    """
    Merge a list of blast output dataframes into a single non-redundant
    dataframe. (See merge_blast_df documentation for details on merge). Grabs
    meta data (isoform, structure, etc.). Converts the "seuqence" column to
    "subject_sequence" and downloads full sequence for that accession from ncbi.

    Parameters
    ----------
    blast_df_list : list
        list of pandas dataframes containing blast output from local_blast,
        ncbi_blast, or read_blast_xml.
    blast_source_list : list, optional
        list of values to put into a "blast_source" column that name where the
        dataframes from blast_df_list came from. Must be the same length as
        blast_df_list. if None, do not populate the "blast_source" column.

    Returns
    -------
    df : pandas.DataFrame
        merged blast dataframe
    """

    # Go through each blast dataframe
    for i in range(len(blast_df_list)):

        # Blast source
        if blast_source_list is not None:
            blast_df_list[i].loc[:,"blast_source"] = blast_source_list[i]

        # Parse the blast output from each line to extract the features useful
        # for downstream analyses -- structure, partial, etc.
        # out_dict will be a dictionary keyed to the new columns we want
        # (structure, etc.) with lists of values as long as the dataframe.
        keep = []
        out_dict = None
        for idx in blast_df_list[i].index:
            parsed = parse_ncbi_line(blast_df_list[i].loc[idx,"title"])

            if parsed is None:
                keep.append(False)
                continue

            if out_dict is None:
                out_dict = {}
                for k in parsed:
                    out_dict[k] = [parsed[k]]
            else:
                for k in parsed:
                    out_dict[k].append(parsed[k])

            keep.append(True)

        # Create mask of goodness
        keep = np.array(keep,dtype=bool)

        # Load newly extracted column into the dataframe
        for k in out_dict:
            blast_df_list[i][k] = pd.NA
            blast_df_list[i].loc[keep,k] = out_dict[k]

        # Drop empty columns
        blast_df_list[i] = blast_df_list[i].loc[keep,:]

    # Drop completely empty blast returns
    blast_df_list = [b for b in blast_df_list if len(b) > 0]
    if len(blast_df_list) == 0:
        err = "BLAST did not return any hits\n"
        raise RuntimeError(err)

    # Merge the dataframes
    df = merge_blast_df(blast_df_list)

    # Download full subject sequences
    downloaded_seq = get_sequences(df.loc[:,"accession"])
    df["subject_sequence"] = df.loc[:,"sequence"]
    df.loc[:,"sequence"] = [s[1] for s in downloaded_seq]

    return df
