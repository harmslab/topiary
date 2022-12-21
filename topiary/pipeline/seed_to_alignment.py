"""
Pipeline that takes a seed dataframe, BLASTS to find sequence hits,
performs quality control, lowers alignment redundancy in a taxonomically
informed fashion, and generates an alignment.
"""

import topiary

from topiary._private import installed
from topiary._private import check
from topiary._private import run_cleanly

import random
import string
import os
import shutil
import urllib.request

def _check_restart(expected_output,restart):
    """
    Decide whether or not to run a given step. 

    Parameters
    ----------
    expected_output : str
        path to the expected output file
    restart : bool
        whether or not we are attempting to restart
    
    Returns
    -------
    run_calc : bool
        whether or not this step needs to be run
    """

    run_calc = True
    if restart:

        # See if json file is there. If so, the run is done.
        if os.path.isfile(expected_output):
            run_calc = False

    return run_calc

def _parse_arguments(out_dir=None,
                     force_species_aware=False,
                     force_not_species_aware=False,
                     restart=False,
                     overwrite=False):
    
    # Make sure the software stack is valid before doing anything
    installed.validate_stack([{"program":"blastp",
                               "min_version":topiary._private.software_requirements["blastp"],
                               "must_pass":True},
                              {"program":"makeblastdb",
                               "min_version":topiary._private.software_requirements["makeblastdb"],
                               "must_pass":True},
                              {"program":"muscle",
                               "min_version":topiary._private.software_requirements["muscle"],
                               "must_pass":True}])

    # --------------------------------------------------------------------------
    # Check sanity of overwrite, restart, and combination

    overwrite = check.check_bool(overwrite,"overwrite")
    restart = check.check_bool(restart,"restart")

    if overwrite and restart:
        err = "overwrite and restart flags are incompatible.\n"
        raise ValueError(err)

    # If no output directory is specified, make up a name
    if out_dir is None:
        if restart:
            err = "To use restart, you must specify an out_dir.\n"
            raise ValueError(err)

        # Make a directory called seed-to-alignment_{dir_counter} where 
        # dir_counter is the first number that does not cause us to overwrite an
        # existing directory
        dir_counter = 1
        out_dir = "seed-to-alignment"
        while os.path.exists(out_dir):
            out_dir = f"seed-to-alignment_{dir_counter:03d}"
            dir_counter += 1

    # --------------------------------------------------------------------------
    # Figure out how to treat species awareness

    force_species_aware = check.check_bool(force_species_aware,"force_species_aware")
    force_not_species_aware = check.check_bool(force_not_species_aware,"force_not_species_aware")
    if force_species_aware and force_not_species_aware:
        err = "force_species_aware and force_not_species_aware cannot both be set to True\n"
        raise ValueError(err)

    # Decide how to treat species awareness
    if force_species_aware:
        species_aware = True
    elif force_not_species_aware:
        species_aware = False
    else:
        species_aware = None

    return overwrite, restart, out_dir, species_aware



@run_cleanly
def seed_to_alignment(seed_df,
                      out_dir=None,
                      seqs_per_column=1,
                      max_seq_number=500,
                      redundancy_cutoff=0.90,
                      worst_align_drop_fx=0.1,
                      sparse_column_cutoff=0.80,
                      align_trim=(0.05,0.95),
                      force_species_aware=False,
                      force_not_species_aware=False,
                      ncbi_blast_db=None,
                      local_blast_db=None,
                      blast_xml=None,
                      move_mrca_up_by=2,
                      local_recip_blast_db=None,
                      min_call_prob=0.95,
                      partition_temp=1,
                      hitlist_size=5000,
                      e_value_cutoff=0.001,
                      gapcosts=(11,1),
                      num_ncbi_blast_threads=1,
                      num_local_blast_threads=-1,
                      restart=False,
                      overwrite=False,
                      keep_recip_blast_xml=False):
    """
    Pipeline that takes a seed dataframe, BLASTs to find sequence hits,
    performs quality control, lowers alignment redundancy in a taxonomically
    informed fashion, and generates an alignment.

    Parameters
    ----------
    seed_df : pandas.DataFrame or str
        Spreadsheet with at least four columns: name, species, sequence, and
        aliases. Can either be a pandas.DataFrame or a string pointing to a 
        spreadsheet file. For details on seed dataframes, see the 
        [documentation](https://topiary-asr.readthedocs.io/en/latest/data_structures.html#seed-dataframe).
    out_dir : str, optional
        Output directory. If not specified, create an output directory with the
        format "seed-to-alignment_{counter}" (where counter increments so 
        previous directories are not overwritten).
    seqs_per_column : float, default=1
        Aim to have this number of sequences per column in the key species
        sequences. (For example, if the key sequence is 100 amino acids long,
        seqs_per_column=1 would aim for 100 sequences; 2 would aim for 200
        sequences).
    max_seq_number : int, default=500
        Maximum number of sequences to get, regardless of seqs_per_column and
        key sequence length.
    redundancy_cutoff : float, default=0.90
        Remove redundant sequences from closely related species is they have
        a sequence identity above redundancy_cutoff.
    worst_align_drop_fx : float, default=0.1
        After alignment, drop approximately this fraction of the sequences,
        selecting those that have long insertions and are missing chunks of
        sequences
    sparse_column_cutoff : float, default=0.80
        When checking alignment quality, a column is sparse if it has gaps in
        more than sparse_column_cutoff sequences.
    align_trim : tuple, default=(0.05,0.95)
        When checking alignment quality, do not score the first and last parts
        of the alignment. Interpreted like a python list slice, but with
        percentages. (0.0,1.0) would not trim; (0.05,0,98) would trim the first
        0.05 off the front and the last 0.02 off the back.

    force_species_aware : bool, default=False
        Lower redundancy in a species-aware fashion, regardless of dataset type.
        By default, topiary will lower sequence redundancy in a species-aware
        fashion for non-microbial datasets, and in a non-species aware fashion
        for microbial datasets. If True, do species aware. 
    force_not_species_aware : bool, default=False
        Lower redundancy in a non-species-aware fashion, regardless of dataset
        type. By default, topiary will lower sequence redundancy in a
        species-aware fashion for non-microbial datasets, and in a non-species
        aware fashion for microbial datasets. If True, do not do species aware. 

    ncbi_blast_db : str, optional
        NCBI blast database to use. (If ncbi_blast_db, local_blast_db and
        blast_xml are all None, ncbi_blast_db is automatically set to "nr").
    local_blast_db : str, optional
        Local blast database to use.
    blast_xml : str or list, optional
        Previously generated blast xml files to load. This argument can be:

         + single xml file (str)
         + list of xml files (list of str)
         + directory (str). Code will grab all .xml files in the directory.

    move_mrca_up_by : int, default=2
        When inferring the phylogenetic context from the seed dataframe, get the
        most recent common ancestor of the seed species, then find the taxonomic
        rank "move_mrca_up_by" levels above that ancestor. For example, if the
        key species all come from marsupials (Theria) and move_mrca_up_by == 2,
        the context will be Amniota (Theria -> Mammalia -> Amniota). Note: if
        all species in the dataset are bacteria or archaea, this argument is
        ignored and the context will be set to Bacteria or Archaea. 
    local_recip_blast_db : str, optional
        Local blast database to use for reciprocal blast. If None, construct a
        reciprocal blast database by downloading the proteomes of the key
        species from the ncbi.
    min_call_prob : float, default=0.95
        Hits from all paralogs that yield a regular expression match to one of
        the aliases from the seed dataframe are weighted by their relative blast
        bit scores. Each paralog is assigned a relative probability. This cutoff
        is the minimum probability the best paralog match must have to result in
        a paralog call. Value should be between 0 and 1 (not inclusive), where
        increasing min_call_prob increases the stringency.
    partition_temp : float, default=1
        When calculating posterior probability of the reciprocal blast paralog
        call, use this for weighting: :code:`2^(bit_score/partition_temp)`.
        :code:`partition_temp` should be a float > 0. A higher value corresponds
        to a higher stringency. (The bit score difference between the best hit
        and the bit scores of other hits would have to be higher to be
        significant). This is a minium value. It may be adjusted automatically
        to avoid numerical problems in the calculation.

    hitlist_size : int, default=5000
        Download only the top :code:`hitlist_size` hits in initial BLAST.
    e_value_cutoff : float, default=0.001
        Only take hits with :code:`e_value` better than :code:`e_value_cutoff`
        in initial blast
    gapcost : tuple, default=(11,1)
        BLAST gapcosts (length 2 tuple of ints) in initial BLAST. The raw score
        of an alignment is the sum of the scores for aligning pairs of residues
        and the scores for gaps. This BLAST search will charge a score
        :code:`gapcost[0]` for the existence of a gap, and the score of
        :code:`gapcost[1]` for each residue in the gap. Thus a gap of :code:`k`
        residues receives a total score of :code:`-(gapcost[0] + gapcost[1]*k)`.
        The default (:code:`(11,1)`) gives scores of :code:`-(11 + 1*k) for each
        gap with :code:`k` residues.
    num_ncbi_blast_threads : int, default=1
        Number of threads to use for NCBI blast. -1 means use all available.
        (The default is 1; multithreading rarely speeds up remote BLAST).
    num_local_blast_threads : int, default=-1
        Number of threads to use for local blast. -1 means all available.

    restart : bool, default=False
        Restart the pipeline from where it stopped in :code:`out_dir`. To use 
        this option, :code:`out_dir` must be specified and point to an existing
        calculation. This option is incompatible with :code:`overwrite = True`.
    overwrite : bool, default=False
        Overwrite :code:`out_dir` if it already exists. This is incompatible
        with :code:`restart = True`.
    keep_recip_blast_xml : bool, default=False
        Whether or not to keep raw blast xml files generated by the pipeline. 

    Returns
    -------
    topiary_dataframe : pandas.DataFrame
        Topiary dataframe with aligned, quality-controlled sequences.
    """

    # Parse the arguments, checking for sanity
    out = _parse_arguments(out_dir=out_dir,
                           force_species_aware=force_species_aware,
                           force_not_species_aware=force_not_species_aware,
                           restart=restart,
                           overwrite=overwrite)

    overwrite, restart, out_dir, species_aware = out

    # --------------------------------------------------------------------------
    # Decide whether to create initial directory

    step_counter = 0
    if issubclass(type(seed_df),str):
        base_name = os.path.split(seed_df)[-1]
    else:
        base_name = "initial-seed-dataframe"
    expected_output = os.path.join(out_dir,f"{step_counter:02d}_{base_name}")
    run_calc = _check_restart(expected_output,restart)

    if run_calc:

        # Try to create the output directory
        if os.path.exists(out_dir):
            if os.path.isdir(out_dir):
                if overwrite:
                    shutil.rmtree(out_dir)
                elif restart:
                    pass
                else:
                    err = f"\nout_dir '{out_dir}' already exists\n\n"
                    raise FileExistsError(err)
            else:
                err = f"\nout_dir '{out_dir}' already exists and is not a directory\n\n"
                raise FileExistsError(err)
        else:
            os.mkdir(out_dir)

        # If the seed_df is a file, copy that into the output directory
        if issubclass(type(seed_df),str):

            # Figure out filename we are going to copy to
            root = f"{step_counter:02d}_{os.path.split(seed_df)[-1]}"
            copy_to = os.path.join(out_dir,root)
    
            # If this is a url, download directly to copy_to
            if seed_df.startswith("https"):
                try:
                    urllib.request.urlretrieve(seed_df,copy_to)
                    seed_df = copy_to
                except Exception as e:
                    err = f"seed_df '{seed_df}' looks like a url, but the file could not be downloaded.\n"
                    raise ValueError(err) from e

            # If the seed_df exists
            if os.path.exists(seed_df):
                if seed_df != copy_to:
                    shutil.copy(seed_df,copy_to)
                seed_df = root

            else:
                err = f"\nseed_df '{seed_df}' not found.\n\n"
                raise FileNotFoundError(err)

    else:
        print(f"Loading existing file {expected_output}.")
        df, key_species, paralog_patterns, species_aware = topiary.io.read_seed(expected_output,species_aware=species_aware)
        seed_df = os.path.split(expected_output)[-1]

    # append paths to local blast resources
    if not blast_xml is None:
        blast_xml = os.path.abspath(blast_xml)

    if not local_blast_db is None:
        local_blast_db = os.path.abspath(local_blast_db)


    # Change into the output directory
    cwd = os.getcwd()
    os.chdir(out_dir)

    step_counter += 1
    expected_output = f"{step_counter:02d}_initial-dataframe.csv"
    run_calc = _check_restart(expected_output,restart)

    if run_calc:

        print("\n-------------------------------------------------------------------")
        print("Building initial topiary dataframe.")
        print("-------------------------------------------------------------------")
        print("",flush=True)

        # If no blast resource specified, default to nr
        if ncbi_blast_db is None and local_blast_db is None and blast_xml is None:
            ncbi_blast_db = "nr"

        kwargs = {"seed_df":seed_df,
                  "ncbi_blast_db":ncbi_blast_db,
                  "local_blast_db":local_blast_db,
                  "blast_xml":blast_xml,
                  "move_mrca_up_by":move_mrca_up_by,
                  "species_aware":species_aware,
                  "hitlist_size":hitlist_size,
                  "e_value_cutoff":e_value_cutoff,
                  "gapcosts":gapcosts,
                  "num_ncbi_blast_threads":num_ncbi_blast_threads,
                  "num_local_blast_threads":num_local_blast_threads,
                  "keep_blast_xml":True}

        out = topiary.df_from_seed(**kwargs)

        df = out[0]
        key_species = out[1]
        paralog_patterns = out[2]
        species_aware = out[3]

        topiary.write_dataframe(df,expected_output)

    else:
        print(f"Loading existing file {expected_output}.")
        df = topiary.read_dataframe(expected_output)

    step_counter += 1
    expected_output = f"{step_counter:02d}_recip-blast-dataframe.csv"
    run_calc = _check_restart(expected_output,restart)

    if run_calc:

        print("\n-------------------------------------------------------------------")
        print("Doing reciprocal blast.")
        print("-------------------------------------------------------------------")
        print("",flush=True)

        # If no reciprocal blast database is specified, construct from key species.
        if local_recip_blast_db is None:

            proteome_list = []
            for k in key_species:
                print(f"Downloading {k} proteome")
                
                # Download proteomes by taxid, not species name because taxid
                # resolver is unambiguous. Species search can lead to no 
                # reference proteome for some species (e.g. E. coli)
                this_taxid = topiary.ncbi.get_taxid(k)
                proteome_list.append(topiary.ncbi.get_proteome(taxid=this_taxid))
                if proteome_list[-1] is None:
                    err = f"\nCould not download proteome for {k} despite multiple\n"
                    err += "attempts. This could be because of high server load\n"
                    err += "at the NCBI or too many requests from your IP address.\n"
                    err += "Try running the pipeline again, appending the --restart\n"
                    err += "flag.\n\n"
                    raise RuntimeError(err)


            base = "".join([random.choice(string.ascii_letters) for _ in range(10)])
            blast_db = f"{base}_local_blast"
            topiary.ncbi.make_blast_db(proteome_list,blast_db,overwrite=True)

            local_recip_blast_db = blast_db

        df = topiary.recip_blast(df,
                                 paralog_patterns=paralog_patterns,
                                 local_blast_db=local_recip_blast_db,
                                 min_call_prob=min_call_prob,
                                 partition_temp=partition_temp,
                                 num_threads=num_local_blast_threads,
                                 keep_blast_xml=keep_recip_blast_xml)

        topiary.write_dataframe(df,expected_output)

    else:
        print(f"Loading existing file {expected_output}.")
        df = topiary.read_dataframe(expected_output)


    # --------------------------------------------------------------------------
    # Redundancy reduction
    
    step_counter += 1
    expected_output = f"{step_counter:02d}_shrunk-dataframe.csv"
    run_calc = _check_restart(expected_output,restart)

    if run_calc:

        print("\n-------------------------------------------------------------------")
        print("Reducing number of sequences.")
        print("-------------------------------------------------------------------")
        print("",flush=True)

        kwargs = {"df":df,
                  "paralog_column":"recip_paralog",
                  "species_tree_aware":species_aware,
                  "seqs_per_column":seqs_per_column*(1 + worst_align_drop_fx),
                  "max_seq_number":max_seq_number*(1 + worst_align_drop_fx),
                  "redundancy_cutoff":redundancy_cutoff,
                  "sparse_column_cutoff":sparse_column_cutoff,
                  "align_trim":align_trim}

        df = topiary.quality.shrink_dataset(**kwargs)
        topiary.write_dataframe(df,expected_output)

    else:
        print(f"Loading existing file {expected_output}.")
        df = topiary.read_dataframe(expected_output)

    step_counter += 1
    expected_output = f"{step_counter:02d}_aligned-dataframe.csv"
    run_calc = _check_restart(expected_output,restart)

    if run_calc:

        print("\n-------------------------------------------------------------------")
        print("Aligning sequences.")
        print("-------------------------------------------------------------------")
        print("",flush=True)

        df = topiary.muscle.align(df)
        topiary.write_dataframe(df,expected_output)
        step_counter += 1

    else:
        print(f"Loading existing file {expected_output}.")

    expected_output = f"{step_counter:02d}_clean-aligned-dataframe.csv"
    run_calc = _check_restart(expected_output,restart)

    if run_calc:

        print("\n-------------------------------------------------------------------")
        print("Polishing alignment and re-aligning.")
        print("-------------------------------------------------------------------")
        print("",flush=True)

        kwargs = {"df":df,
                  "fx_sparse_percentile":(1 - worst_align_drop_fx),
                  "sparse_run_percentile":(1 - worst_align_drop_fx),
                  "fx_missing_percentile":(1 - worst_align_drop_fx),
                  "realign":True,
                  "sparse_column_cutoff":sparse_column_cutoff}

        df = topiary.quality.polish_alignment(**kwargs)
        topiary.write_dataframe(df,expected_output)
        pretty_name = os.path.join(out_dir,expected_output)
        step_counter += 1

        topiary.write_fasta(df,
                            f"{step_counter:02d}_alignment.fasta",
                            seq_column="alignment",
                            label_columns=["species","recip_paralog"])

        print("\n-------------------------------------------------------------------")
        print(f"Dataset in {pretty_name}.")
        print("-------------------------------------------------------------------")
        print("",flush=True)

    else:
        print(f"Loading existing file {expected_output}.")

    os.chdir(cwd)

    return df
