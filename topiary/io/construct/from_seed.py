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
                 blast_xml=None,
                 move_mrca_up_by=2,
                 hitlist_size=5000,
                 e_value_cutoff=0.001,
                 gapcosts=(11,1),
                 num_threads=1,
                 keep_blast_xml=False,
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
    blast_xml : str or list, optional
        previously generated blast xml files to load. This argument can be:

         + single xml file (str)
         + list of xml files (list of str)
         + directory (str). Code will grab all .xml files in the directory.

    move_mrca_up_by : int, default=2
        when inferring the phylogenetic context from the seed dataframe, get the
        most recent common ancestor of the seed species, then find the taxonomic
        rank "move_mrca_up_by" levels above that ancestor. For example, if the
        key species all come from marsupials (Theria) and move_mrca_up_by == 2,
        the context will be Amniota (Theria -> Mammalia -> Amniota).
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
    keep_blast_xml : bool, default=False
        whether or not to keep raw blast xml output
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
    key_species : numpy.array
        list if key species to keep during the analysis
    paralog_patterns : list
        list of compiled regular expressions (extracted from aliases) to use to
        try to match paralogs.

    Notes
    -----
    Every sequence in the original seed dataframe will have :code:`always_keep`
    set to :code:`True`, so they will not be deleted by subsequent quality
    control steps.
    """

    # Load the seed dataframe
    seed_df = topiary.io.load_seed_dataframe(seed_df)

    # -----------------------------------------------------------------------
    # Get key_species and paralog_patterns dict

    key_species = np.unique(seed_df.loc[seed_df.loc[:,"key_species"],"species"])
    key_species.sort()

    paralog_patterns = {}
    for i, k in enumerate(seed_df.loc[:,"name"]):
        paralog_patterns[k] = seed_df.loc[seed_df.index[i],"aliases"]

    # Validate blast database arguments
    if ncbi_blast_db is not None and local_blast_db is not None and blast_xml is None:
        err = "\nYou must specify ncbi_blast_db or local_blast_db or\n"
        err += "blast_xml or some combination thereof to find sequences.\n\n"
        raise ValueError(err)

    # database, hitlist_size, e_value_cutoff, gapcosts, num_threads, and
    # kwargs will all be validated by the blast call itself.

    # ncbi blast
    blast_df = []
    if ncbi_blast_db is not None:

        # Infer phylogenetic context from key species
        phylo_context = topiary.opentree.get_mrca(key_species,
                                                  move_up_by=move_mrca_up_by)
        try:
            taxid = phylo_context["taxid"]
        except KeyError:
            taxid = None

        tmp_blast_df = topiary.ncbi.ncbi_blast(seed_df.sequence,
                                               db="nr",
                                               taxid=taxid,
                                               hitlist_size=hitlist_size,
                                               e_value_cutoff=e_value_cutoff,
                                               gapcosts=gapcosts,
                                               num_threads=num_threads,
                                               keep_blast_xml=keep_blast_xml,
                                               **kwargs)
        blast_df.extend(tmp_blast_df)

    # local blast
    if local_blast_db is not None:

        tmp_blast_df = topiary.ncbi.local_blast(seed_df.sequence,
                                                db=local_blast_db,
                                                hitlist_size=hitlist_size,
                                                e_value_cutoff=e_value_cutoff,
                                                gapcosts=gapcosts,
                                                num_threads=num_threads,
                                                keep_blast_xml=keep_blast_xml,
                                                **kwargs)
        blast_df.extend(tmp_blast_df)

    # Load blast xml
    if blast_xml is not None:
        err = "Not implemented\n"
        raise NotImplementedError(err)


    print("Parsing BLAST output",flush=True)

    # blast_df is a list of dataframes, one for each query in seed_df.sequence
    for i, d in enumerate(blast_df):
        d.to_csv(f"blast_{i}.csv")

    # Go through each blast dataframe
    for i in range(len(blast_df)):

        if len(blast_df[i]) == 0:
            seed_name = seed_df.loc[seed_df.index[i],"name"]
            seed_species = seed_df.loc[seed_df.index[i],"species"]
            print(f"There were no BLAST hits for {seed_name} from {seed_species}",
                  flush=True)
            continue

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

    # Will convert to topiary dataframe
    df = check.check_topiary_dataframe(df)
    df = pd.concat((seed_df,df),ignore_index=True)

    # Set always_keep and key_species for new hits
    df.loc[pd.isna(df["always_keep"]),"always_keep"] = False
    df.loc[pd.isna(df["key_species"]),"key_species"] = False

    print("Getting OTT species ids for all species.",flush=True)

    # Get ott id for all sequences, setting False for those that can't be
    # found/resolved
    df = topiary.get_ott(df,phylo_context=phylo_context,verbose=False)

    # Create nicknames for sequences in dataframe
    df = topiary.create_nicknames(df,paralog_patterns)

    return df, phylo_context, key_species, paralog_patterns
