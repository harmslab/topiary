"""
Pipeline that takes a seed dataframe, BLASTS to find sequence hits,
performs quality control, lowers alignment redundancy in a taxonomically
informed fashion, and generates an alignment.
"""

import topiary

from topiary._private import installed
from topiary._private import check

import numpy as np
import pandas as pd

import random, string, os, shutil

def _check_restart(expected_output,restart):

    run_calc = True
    if restart:

        # See if json file is there. If so, the run is done.
        if os.path.isfile(expected_output):
            run_calc = False

    return run_calc

def seed_to_alignment(seed_df,
                      out_dir,
                      seqs_per_column=1,
                      max_seq_number=500,
                      redundancy_cutoff=0.90,
                      worst_align_drop_fx=0.1,
                      sparse_column_cutoff=0.80,
                      align_trim=(0.05,0.95),
                      ncbi_blast_db="nr",
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
                      keep_recip_blast_xml=False,
                      verbose=False):
    """
    Pipeline that takes a seed dataframe, BLASTs to find sequence hits,
    performs quality control, lowers alignment redundancy in a taxonomically
    informed fashion, and generates an alignment.

    Parameters
    ----------
    seed_df : pandas.DataFrame or str
        dataframe with at least four columns: name, species, sequence,
        and aliases. See documentation on seed dataframes for details.
    out_dir : str
        output directory

    seqs_per_column : float, default=1
        aim to have this number of sequences per column in the key species
        sequences. (For example, if the key sequence is 100 amino acids long,
        seqs_per_column=1 would aim for 100 sequences; 2 would aim for 200
        sequences).
    max_seq_number : int, default=500
        maximum number of sequences to get, regardless of seqs_per_column and
        key sequence length.
    redundancy_cutoff : float, default=0.90
        merge sequences from closely related species with sequence identity
        above cutoff.
    worst_align_drop_fx : float, default=0.1
        after alignment, drop approximately this fraction of the sequences,
        selecting those that have long insertions and are missing chunks of
        sequences
    sparse_column_cutoff : float, default=0.80
        when checking alignment quality, a column is sparse if it has gaps in
        more than sparse_column_cutoff sequences.
    align_trim : tuple, default=(0.05,0.95)
        when checking alignment quality, do not score the first and last parts
        of the alignment. Interpreted like a slice, but with percentages.
        (0.0,1.0) would not trim; (0.05,0,98) would trim the first 0.05 off the
        front and the last 0.02 off the back.

    ncbi_blast_db : str or None, default="nr"
        NCBI blast database to use.
    local_blast_db : str, optional
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
        the context will be Amniota (Theria -> Mammalia -> Amniota).
    local_recip_blast_db : str, optional
        Local blast database to use for reciprocal blast. If None, construct a
        reciprocal blast database by downloading the proteomes of the key
        species from the ncbi.
    min_call_prob : float, default=0.95
        hits from all paralogs that yield a regular expression match to one of
        the aliases from the seed dataframe are weighted by their relative blast
        bit scores. Each paralog is assigned a relative probability. This cutoff
        is the minimum probability the best paralog match must have to result in
        a paralog call. Value should be between 0 and 1 (not inclusive), where
        increasing min_call_prob increases the stringency.
    partition_temp : float, default=1
        when calculating posterior probability of the reciprocal blast paralog
        call, use this for weighting: 2^(bit_score/partition_temp).
        partition_temp should be a float > 0. A higher value corresponds to a
        higher stringency. (The bit score difference between the best hit and
        the bit scores of other hits would have to be higher to be significant).
        This is a minium value. It may be adjusted automatically to avoid
        numerical problems in the calculation.

    hitlist_size : int, default=5000
        download only the top hitlist_size hits in initial blast
    e_value_cutoff : float, default=0.001
        only take hits with e_value better than e_value_cutoff in initial blast
    gapcost : tuple, default=(11,1)
        BLAST gapcosts (length 2 tuple of ints) in initial blast
    num_ncbi_blast_threads : int, default=1
        number of threads to use for NCBI blast. -1 means use all available.
        (Multithreading rarely speeds up remote BLAST).
    num_local_blast_threads : int, default=-1
        number of threads to use for local blast. -1 means all available.

    restart : bool, default=False
        restart job from where it stopped in output directory. incompatible with
        overwrite
    overwrite : bool, default=False
        overwrite out_dir if it already exists. incompatible with restart
    keep_recip_blast_xml : bool, default=False
        whether or not to keep raw blast xml output
    verbose : bool, default=False
        verbosity of output

    Returns
    -------
    topiary_dataframe : pandas.DataFrame
        Topiary dataframe with aligned, quality-controlled sequences.
    """

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

    out_dir = str(out_dir)

    step_counter = 0
    expected_output = os.path.join(out_dir,f"{step_counter:02d}_{os.path.split(seed_df)[-1]}")
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
        if type(seed_df) is str:
            if os.path.exists(seed_df):
                root = f"{step_counter:02d}_{os.path.split(seed_df)[-1]}"
                copy_to = os.path.join(out_dir,root)
                shutil.copy(seed_df,copy_to)
                seed_df = root
            else:
                err = f"\nseed_df '{seed_df}' not found.\n\n"
                raise FileNotFoundError(err)

    else:
        print(f"Loading existing file {expected_output}.")
        df, key_species, paralog_patterns = topiary.io.read_seed(expected_output)
        seed_df = os.path.split(expected_output)[-1]

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

        kwargs = {"seed_df":seed_df,
                  "ncbi_blast_db":ncbi_blast_db,
                  "local_blast_db":local_blast_db,
                  "blast_xml":blast_xml,
                  "move_mrca_up_by":move_mrca_up_by,
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
                proteome_list.append(topiary.ncbi.get_proteome(species=k))
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
                                 num_threads=num_local_blast_threads,
                                 keep_blast_xml=keep_recip_blast_xml)

        topiary.write_dataframe(df,expected_output)

    else:
        print(f"Loading existing file {expected_output}.")
        df = topiary.read_dataframe(expected_output)

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
                  "fx_sparse_percential":(1 - worst_align_drop_fx),
                  "sparse_run_percentile":(1 - worst_align_drop_fx),
                  "fx_missing_percentile":(1 - worst_align_drop_fx),
                  "realign":True,
                  "sparse_column_cutoff":sparse_column_cutoff}

        df = topiary.quality.polish_alignment(**kwargs)
        topiary.write_dataframe(df,expected_output)
        step_counter += 1

        topiary.write_fasta(df,
                            f"{step_counter:02d}_alignment.fasta",
                            seq_column="alignment",
                            label_columns=["species","recip_paralog"])

    else:
        print(f"Loading existing file {expected_output}.")

    os.chdir(cwd)

    return df
