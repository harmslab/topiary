"""
Construct a topiary dataframe from a seed dataframe.
"""

import topiary
from topiary._private import check
from topiary.external.ncbi.blast import merge_and_annotate

import numpy as np
import pandas as pd

import re

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

    Notes
    -----
    Every sequence in the original seed dataframe will have :code:`always_keep`
    set to :code:`True`, so they will not be deleted by subsequent quality
    control steps.
    """

    # Load the seed dataframe
    seed_df, key_species, paralog_patterns = topiary.io.read_seed(seed_df)

    # Validate proteome
    no_proteome_key = []
    for k in key_species:
        proteome_id, err = topiary.ncbi.entrez.get_proteome_ids(species=k)
        if proteome_id is None:
            no_proteome_key.append((k,err))

    if len(no_proteome_key) > 0:
        err = "\nCould not download proteomes from the NCBI assemblies database\n"
        err += "for all key species in the seed dataframe. Proteomes are required\n"
        err += "for reciprocal BLAST. Please check the spelling of your species\n"
        err += "and/or select species with avaiable proteomes. The species that\n"
        err += "raised errors are:\n\n"
        for n in no_proteome_key:
            err += f"{n[0]} ({n[1].strip()})\n"
        err += "\n"

        raise RuntimeError(err)

    # Validate blast database arguments
    if ncbi_blast_db is not None and local_blast_db is not None and blast_xml is None:
        err = "\nYou must specify ncbi_blast_db or local_blast_db or\n"
        err += "blast_xml or some combination thereof to find sequences.\n\n"
        raise ValueError(err)

    # database, hitlist_size, e_value_cutoff, gapcosts, num_threads, and
    # kwargs will all be validated by the blast call itself.

    # ncbi blast
    blast_df = []
    blast_source = []
    if ncbi_blast_db is not None:

        # Infer phylogenetic context from key species
        phylo_context = topiary.opentree.ott_mrca(species_list=key_species,
                                                  move_up_by=move_mrca_up_by)
        try:
            taxid = phylo_context["taxid"]
        except KeyError:
            taxid = None

        tmp_blast_df = topiary.ncbi.ncbi_blast(seed_df.sequence,
                                               db=ncbi_blast_db,
                                               taxid=taxid,
                                               hitlist_size=hitlist_size,
                                               e_value_cutoff=e_value_cutoff,
                                               gapcosts=gapcosts,
                                               num_threads=num_threads,
                                               keep_blast_xml=keep_blast_xml,
                                               **kwargs)

        for i in range(len(tmp_blast_df)):
            if len(tmp_blast_df[i]) == 0:
                seed_name = seed_df.loc[seed_df.index[i],"name"]
                seed_species = seed_df.loc[seed_df.index[i],"species"]

                w = f"\nThere were no hits from the ncbi {ncbi_blast_db} database\n"
                w += f"using {seed_name} from {seed_species}.\n\n"
                print(w,flush=True)

        blast_df.extend(tmp_blast_df)
        for idx in seed_df.index:
            db = ncbi_blast_db
            name = seed_df.loc[idx,"name"]
            species = seed_df.loc[idx,"species"]
            blast_source.append(f"ncbi {db}|{name}|{species}")

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

        for i in range(len(tmp_blast_df)):
            if len(tmp_blast_df[i]) == 0:
                seed_name = seed_df.loc[seed_df.index[i],"name"]
                seed_species = seed_df.loc[seed_df.index[i],"species"]

                w = f"\nThere were no hits from the local {local_blast_db} database\n"
                w += f"using {seed_name} from {seed_species}.\n\n"
                print(w,flush=True)

        blast_df.extend(tmp_blast_df)
        for idx in seed_df.index:
            db = local_blast_db
            name = seed_df.loc[idx,"name"]
            species = seed_df.loc[idx,"species"]
            blast_source.append(f"local {db}|{name}|{species}")

    # Load blast xml
    if blast_xml is not None:

        tmp_blast_df, xml_files = topiary.io.read_blast_xml(blast_xml)

        blast_df.extend(tmp_blast_df)
        for x in xml_files:
            blast_source.append(f"xml {x}")

    # Merge and annotate the blast input from all sources
    df = merge_and_annotate(blast_df,blast_source)

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
    df = topiary.get_ott(df,verbose=False)

    # Create nicknames for sequences in dataframe
    df = topiary.create_nicknames(df,paralog_patterns)

    return df, key_species, paralog_patterns
