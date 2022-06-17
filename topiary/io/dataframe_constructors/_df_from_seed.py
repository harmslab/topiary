__description__ = \
"""
Construct a topiary dataframe from a seed spreadsheet.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-09"

import topiary
from topiary import _arg_processors

import numpy as np
import pandas as pd

def df_from_seed(seed_df,
                 ncbi_blast_db="nr",
                 local_blast_db=None,
                 hitlist_size=5000,
                 e_value_cutoff=0.001,
                 gapcosts=(11,1),
                 num_threads=1,
                 **kwargs):


    out = topiary.io.load_seed_dataframe(seed_df)
    seed_df = out[0]
    phylo_context = out[1]
    key_species = out[2]
    paralog_patterns = out[3]

    # Validate dataframe
    seed_df = _arg_processors.process_topiary_dataframe(seed_df)

    # Validate blast database arguments
    if ncbi_blast_db is None and local_blast_db is None:
        err = "\nPlease specificy either ncbi_blast_db OR local_blast_db\n\n"
        raise ValueError(err)

    if ncbi_blast_db is not None and local_blast_db is not None:
        err = "\nPlease specificy either ncbi_blast_db OR\n"
        err += "local_blast_db, but not both.\n\n"
        raise ValueError(err)

    # database, hitlist_size, e_value_cutoff, gapcosts, num_threads, and
    # kwargs will all be validated by the blast call itself.

    # ncbi blast
    if ncbi_blast_db:

        # Get the taxid to limit the search to if phylo_context is given
        if phylo_context is not None:
            taxid = topiary.opentree.phylo_to_taxid[phylo_context]
        else:
            taxid = None

        blast_df = topiary.ncbi.ncbi_blast(seed_df.sequence,
                                           db="nr",
                                           taxid=taxid,
                                           hitlist_size=hitlist_size,
                                           e_value_cutoff=e_value_cutoff,
                                           gapcosts=gapcosts,
                                           num_threads=num_threads,
                                           **kwargs)

    # local blast
    if local_blast_db:

        blast_df = topiary.ncbi.local_blast(seed_df.sequence,
                                            db=local_blast_db,
                                            hitlist_size=hitlist_size,
                                            e_value_cutoff=e_value_cutoff,
                                            gapcosts=gapcosts,
                                            num_threads=num_threads,
                                            **kwargs)

    # blast_df is a list of dataframes, one for each hit in seed_df.sequence

    # Go through each blast dataframe
    for i, this_df in enumerate(blast_df):

        # Assign the blast query to a useful name (i.e. LY96|Homo sapiens)
        this_df.loc[:,"query"] = f"{seed_df.name.iloc[i]}|{seed_df.species.iloc[i]}"

        # Parse the blast output from each line to extract the features useful
        # for downstream analyses -- structure, partial, etc.
        # out_dict will be a dictionary keyed to the new columns we want
        # (structure, etc.) with lists of values as long as the dataframe.
        out_dict = None
        for idx in this_df.index:
            parsed = topiary.external.ncbi.parse_ncbi_line(this_df.loc[idx,"title"])
            if out_dict is None:
                out_dict = {}
                for k in parsed:
                    out_dict[k] = [parsed[k]]
            else:
                for k in parsed:
                    out_dict[k].append(parsed[k])

        # Load newly extracted column into the dataframe
        for k in out_dict:
            this_df[k] = out_dict[k]


    df = topiary.ncbi.merge_blast_df(blast_df)

    # Append uid (to prevent user-visible topiary warnings) and then
    # convert to a topiary dataframe
    df["uid"] = topiary._private.generate_uid(len(df.sequence))

    df = topiary._arg_processors.process_topiary_dataframe(df)

    df = pd.concat((seed_df,df),ignore_index=True)

    # Set new hits to always_keep
    df.loc[pd.isna(df["always_keep"]),"always_keep"] = False

    # Get ott id for all sequences, setting False for those that can't be
    # found/resolved
    df = topiary.get_ott_id(df,phylo_context=phylo_context,verbose=False)

    # Create nicknames for sequences in dataframe
    df = topiary.create_nicknames(df,paralog_patterns)

    return df, phylo_context, key_species, paralog_patterns
