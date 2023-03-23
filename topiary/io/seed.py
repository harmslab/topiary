"""
Functions for working with seed dataframes
"""

import topiary
from topiary._private import check
from topiary.opentree import species_to_ott
from topiary.opentree import ott_to_resolvable
from topiary.ncbi.blast import merge_and_annotate
from topiary.ncbi.blast import read_blast_xml
from topiary.io.paralog_patterns import load_paralog_patterns

import numpy as np
import pandas as pd

import warnings

def read_seed(df,
              species_aware=None):
    """
    Read a seed data frame and extract alias patterns and key species.

    Parameters
    ----------
    df: pandas.DataFrame or str
        seed dataframe containing seed sequences to launch the analysis. df can
        be a pandas dataframe or a string pointing to a spreadsheet file.
    species_aware : bool or None, default=None
        Whether or not read seed in a species-aware fashion. If True, require
        all species be resolvable in the seed dataset. If False, do not require
        resolvable. If None, choose automatically. (If microbial, set to False;
        if not microbial, set to True). 

    Returns
    -------
    topiary_dataframe : pandas.DataFrame
        new topiary dataframe built from the seed dataframe
    key_species : numpy.array
        list if key species to keep during the analysis
    paralog_patterns : list
        list of compiled regular expressions to use to try to match paralogs
    species_aware : bool
        whether or not this should be treated as a species aware calculation

    Notes
    -----

    The seed dataframe is expected to have at least four columns:

    + :code:`species`: species names for seed sequences in binomial format (i.e.
      Homo sapiens or Mus musculus)
    + :code:`name`: name of each sequence (i.e. LY96)
    + :code:`aliases`: other names for this sequence found in different
      databases/species, separated by ; (i.e. LY96;MD2;ESOP1)
    + :code:`sequence`: amino acid sequences for these proteins.

    It may have one other optional column:

    + :code:`key_species`: True/False. Indicates whether or not this species
      should be used as a key species for reciprocal BLASTing.

    Other columns in the dataframe are kept but not used by topiary.
    """

    # -----------------------------------------------------------------------
    # Load dataframe and check it's sanity

    if type(df) is not type(pd.DataFrame()):

        if issubclass(type(df),str):
            extension = df.split(".")[-1].strip().lower()
            if extension in ["xls","xlsx"]:
                df = pd.read_excel(df)
            elif extension == "csv":
                df = pd.read_csv(df,sep=",")
            elif extension == "tsv":
                df = pd.read_csv(df,sep="\t")
            else:
                df = pd.read_csv(df,sep=None,engine="python")
        else:
            err = f"Could not figure out how to read df '{df}'\n\n"
            raise ValueError(err)

    required_columns = ["species","name","aliases","sequence"]
    df.columns = [c.lower().strip() for c in df.columns]
    for c in required_columns:

        try:
            df[c]
        except KeyError:
            err = f"\nDataframe must have a {c} column. Required columns are:\n"
            for r in required_columns:
                err += f"    {r}\n"
            err += "\n"
            raise ValueError(err)

        # Strip extra leading/trailing white space
        df.loc[:,c] = df.loc[:,c].str.strip()

    # -----------------------------------------------------------------------
    # Get species ott

    bad_species = []
    ott_list, species_list, _ = species_to_ott(np.unique(df.loc[:,"species"]))
    for i in range(len(species_list)):
        if ott_list[i] is None:
            bad_species.append(species_list[i])

    # -------------------------------------------------------------------------
    # Figure out whether to keep_unresolvable if not specified in function call

    if species_aware is None:

        if len(bad_species) > 0:
            w = "Warning. Topiary could not find all species in the input dataframe\n"
            w += "on the Open Tree of Life database. Topiary will determine whether\n"
            w += "this is a purely microbial protein or not based on the species\n"
            w += "that could be identified.\n"
            w += "This unrecognized species are:\n\n"
            for b in bad_species:
                w += f"    '{b}'\n"
            warnings.warn(w)


        tmp_ott_list = [ott for ott in ott_list if ott is not None]

        # Figure out if the default is to reconcile or not based on taxonomic 
        # distribution of alignment. 
        mrca = topiary.opentree.ott_to_mrca(ott_list=tmp_ott_list,
                                            avoid_all_life=True)

        if mrca["is_microbial"]:
            species_aware = False
        else:
            species_aware = True

    # -------------------------------------------------------------------------
    # If we must resolve all species, do some sanity checks

    if species_aware:

        if len(bad_species) > 0:
            err = "\nNot all input species were found in the Open Tree of Life\n"
            err += "database. To troubleshoot the problem, you can visit\n"
            err += "https://tree.opentreeoflife.org/taxonomy/browse and search for\n"
            err += "the following species manually:\n\n"
            for b in bad_species:
                err += f"    '{b}'\n"
            raise ValueError(err)

        # Make sure all input species can be resolved on the OTT synthetic tree
        resolved = ott_to_resolvable(ott_list)
        unresolvable_species = []
        for i in range(len(resolved)):
            if not resolved[i]:
                unresolvable_species.append((species_list[i],ott_list[i]))

        if len(unresolvable_species) > 0:
            err = "\nNot all input species can be resolved on the latest Open Tree of\n"
            err += "Life synthetic tree. To troubleshoot the problem, you can visit\n"
            err += "https://tree.opentreeoflife.org/taxonomy/browse and search for\n"
            err += "the OTT accession numbers of these species:\n"

            for b in unresolvable_species:
                err += f"    '{b[0]}' (OTT '{b[1]}')\n"
            raise ValueError(err)

    # -----------------------------------------------------------------------
    # Get key_species

    # Look for key_species column. If present, make sure it is bool. If not, add
    # it and set all to True
    if "key_species" in df.columns:
        df.loc[:,"key_species"] = check.column_to_bool(df.loc[:,"key_species"],
                                                       "key_species")
    else:
        df["key_species"] = True

    key_species = np.unique(df.loc[df.loc[:,"key_species"],"species"])
    key_species.sort()

    # -----------------------------------------------------------------------
    # Get regular expression for searching for paralogs

    alias_dict = {}
    for idx in df.index:

        # Construct list of all aliases given by the user for a given name
        key = str(df.loc[idx,"name"]).strip()
        values = [v.strip().lower() for v in df.loc[idx,"aliases"].split(";")]

        # Get rid of ""
        values = [v for v in values if v != ""]

        try:
            alias_dict[key].extend(values)
        except KeyError:
            alias_dict[key] = values[:]

    paralog_patterns = load_paralog_patterns(alias_dict)

    # Load newly created paralog patterns into aliases column
    new_aliases = []
    for idx in df.index:

        name = str(df.loc[idx,"name"]).strip()
        new_aliases.append(paralog_patterns[name].pattern)

    df.loc[:,"aliases_regex"] = new_aliases


    # -----------------------------------------------------------------------
    # Convert dataframe into a topiary dataframe

    # add always keep column to this dataframe
    df["always_keep"] = True

    # Drop on uid column to avoid warning when topiary creates in arg_processor
    try:
        df["uid"]
    except KeyError:
        df["uid"] = topiary._private.generate_uid(len(df["always_keep"]))

    df = check.check_topiary_dataframe(df)

    return df, key_species, paralog_patterns, species_aware


def df_from_seed(seed_df,
                 ncbi_blast_db="nr",
                 local_blast_db=None,
                 blast_xml=None,
                 move_mrca_up_by=2,
                 species_aware=None,
                 hitlist_size=5000,
                 e_value_cutoff=0.001,
                 gapcosts=(11,1),
                 num_ncbi_blast_threads=1,
                 num_local_blast_threads=-1,
                 keep_blast_xml=False,
                 **kwargs):
    """
    Construct a topiary dataframe from a seed dataframe, blasting to fill in the
    sequences. This can blast an NCBI database, local database, and/or read in
    previously-run blast xml files.

    Parameters
    ----------
    seed_df : pandas.DataFrame or str
        seed dataframe containing seed sequences to launch the analysis. df can
        be a pandas dataframe or a string pointing to a spreadsheet file.
    ncbi_blast_db : str or None, default="nr"
        NCBI blast database to use.
    local_blast_db : str or None, default=None
        Local blast database to use.
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
        the context will be Amniota (Theria -> Mammalia -> Amniota). Note: If 
        the seed dataframe consists entirely of Bacterial or Archaeal sequences,
        the mrca will be set to the appropriate domain, not a local species
        ancestors.
    species_aware : bool or None, optional
        If True, do analysis in species-aware fashion; if False, ignore species; 
        if None, infer this from the dataset. (Microbial datasets will be False;
        non-microbial datasets will be True.)

    hitlist_size : int, default=5000
        download only the top hitlist_size hits
    e_value_cutoff : float, default=0.001
        only take hits with e_value better than e_value_cutoff
    gapcost : tuple, default=(11,1)
        BLAST gapcosts (length 2 tuple of ints)
    num_ncbi_blast_threads : int, default=1
        number of threads to use for NCBI blast. -1 means use all available.
        (Multithreading rarely speeds up remote BLAST).
    num_local_blast_threads : int, default=-1
        number of threads to use for local blast. -1 means all available.
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
    key_species : numpy.array
        list if key species to keep during the analysis
    paralog_patterns : list
        list of compiled regular expressions to use to try to match paralogs
    species_aware : bool
        whether or not this should be treated as a species aware calculation

    Notes
    -----
    Every sequence in the original seed dataframe will have :code:`always_keep`
    set to :code:`True`, so they will not be deleted by subsequent quality
    control steps.
    """

    # Load the seed dataframe
    seed_df, key_species, paralog_patterns, species_aware = topiary.io.read_seed(seed_df,
                                                                                 species_aware=species_aware)

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

        print(f"BLASTing against NCBI database {ncbi_blast_db}")

        # Infer phylogenetic context from key species
        phylo_context = topiary.opentree.ott_to_mrca(species_list=key_species,
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
                                               num_threads=num_ncbi_blast_threads,
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

        print(f"BLASTing against local database {local_blast_db}")

        tmp_blast_df = topiary.ncbi.local_blast(seed_df.sequence,
                                                db=local_blast_db,
                                                hitlist_size=hitlist_size,
                                                e_value_cutoff=e_value_cutoff,
                                                gapcosts=gapcosts,
                                                num_threads=num_local_blast_threads,
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

        print(f"Loading existing blast results from from {blast_xml}")

        tmp_blast_df, xml_files = read_blast_xml(blast_xml)

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

    # Drop any "synthetic" sequences that came in. (These actually have an OTT
    # and are placed as an outgroup to all life!)
    synth_mask = df.species.str.match("synthetic")
    df.loc[synth_mask,"keep"] = False

    # Combine seed and downloaded sequences.
    df = pd.concat((seed_df,df),ignore_index=True)

    # Set always_keep and key_species for new hits
    df.loc[pd.isna(df["always_keep"]),"always_keep"] = False
    df.loc[pd.isna(df["key_species"]),"key_species"] = False

    print("Getting OTT species ids for all species.",flush=True)

    # Get ott id for all sequences, setting False for those that can't be
    # found/resolved (unless we're not species aware)
    if species_aware: 
        keep_anyway = False
    else:
        keep_anyway = True
    df = topiary.get_df_ott(df,verbose=False,keep_anyway=keep_anyway)

    # Create nicknames for sequences in dataframe
    df = topiary.create_nicknames(df,paralog_patterns)

    return df, key_species, paralog_patterns, species_aware
