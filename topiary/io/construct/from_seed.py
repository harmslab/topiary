"""
Construct a topiary dataframe from a seed dataframe.
"""

import topiary
from topiary._private import check

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
    """
    Construct a topiary dataframe from a seed dataframe, blasting to fill in the
    sequences.

    Parameters
    ----------
    seed_df : pandas.DataFrame or str
        seed dataframe containing seed sequences to launch the analysis. df can
        be a pandas dataframe or a string pointing to a spreadsheet file.
    ncbi_blast_db : str or None, default="nr"
        NCBI blast database to use. If None, use a local database. Incompatible
        with local_blast_db.
    local_blast_db : str or None, default=None
        Local blast database to use. If None, use an NCBI database. Incompatible
        with ncbi_blast_db.
    hitlist_size : int, default=5000
        download only the top hitlist_size hits
    e_value_cutoff : float, default=0.001
        only take hits with e_value better than e_value_cutoff
    gapcost : tuple, default=(11,1)
        BLAST gapcosts (length 2 tuple of ints)
    num_threads : int, default=1
        number of threads to use for BLAST. If -1, use all available threads.
        Note: multithreading rarely speeds up NCBI blast queries, but can
        dramatically speed up local blast searches.
    **kwargs : dict, optional
        extra keyword arguments are passed directly to biopython
        NcbiblastXXXCommandline (for local blast) or qblast (for remote
        blast). These take precedence over anything specified above
        (hitlist_size, for example).

    Returns
    -------
    topiary_dataframe : pandas.DataFrame
        topiary dataframe with sequences found from seed sequence.
    phylo_context : str
        string descriptor of the phylogenetic context on opentreeoflife.org
        (i.e. "Animals" or "All life")
    key_species : list
        list if key species to keep during the analysis
    paralog_patterns : list
        list of compiled regular expressions (extracted from aliases) to use to
        try to match paralogs.

    Notes
    -----
    Every sequence in the original seed dataframe will have `always_keep` set
    to True, so they will not be deleted by subsequent quality control steps.
    """

    out = topiary.io.load_seed_dataframe(seed_df)
    seed_df = out[0]
    phylo_context = out[1]
    key_species = out[2]
    paralog_patterns = out[3]

    # Validate dataframe
    seed_df = check.check_topiary_dataframe(seed_df)

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
    for i in range(len(blast_df)):

        # Assign the blast query to a useful name (i.e. LY96|Homo sapiens)
        blast_df[i].loc[:,"query"] = f"{seed_df.name.iloc[i]}|{seed_df.species.iloc[i]}"

        # Parse the blast output from each line to extract the features useful
        # for downstream analyses -- structure, partial, etc.
        # out_dict will be a dictionary keyed to the new columns we want
        # (structure, etc.) with lists of values as long as the dataframe.
        keep = []
        out_dict = None
        for idx in blast_df[i].index:
            parsed = topiary.external.ncbi.parse_ncbi_line(blast_df[i].loc[idx,"title"])

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
            blast_df[i][k] = pd.NA
            blast_df[i].loc[keep,k] = out_dict[k]

        # Drop empty columns
        blast_df[i] = blast_df[i].loc[keep,:]

    # Drop completely empty blast returns
    blast_df = [b for b in blast_df if len(b) > 0]
    if len(blast_df) == 0:
        err = "BLAST did not return any hits\n"
        raise RuntimeError(err)

    df = topiary.ncbi.merge_blast_df(blast_df)

    # Append uid (to prevent user-visible topiary warnings) and then
    # convert to a topiary dataframe
    df["uid"] = topiary._private.generate_uid(len(df.sequence))

    df = check.check_topiary_dataframe(df)

    df = pd.concat((seed_df,df),ignore_index=True)

    # Set new hits to always_keep
    df.loc[pd.isna(df["always_keep"]),"always_keep"] = False

    # Get ott id for all sequences, setting False for those that can't be
    # found/resolved
    df = topiary.get_ott_id(df,phylo_context=phylo_context,verbose=False)

    # Create nicknames for sequences in dataframe
    df = topiary.create_nicknames(df,paralog_patterns)

    return df, phylo_context, key_species, paralog_patterns
