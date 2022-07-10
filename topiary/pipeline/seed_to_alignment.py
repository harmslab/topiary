"""
Pipeline that takes a seed dataframe, BLASTS to find sequence hits,
performs quality control, lowers alignment redundancy in a taxonomically
informed fashion, and generates an alignment.
"""

import topiary

from topiary._private import installed

import numpy as np
import pandas as pd

import random, string, os, shutil

def seed_to_alignment(seed_df,
                      out_dir,
                      target_seq_number=500,
                      ncbi_blast_db="nr",
                      local_blast_db=None,
                      local_recip_blast_db=None,
                      within_species_redundancy_cutoff=0.99,
                      sparse_column_cutoff=0.95,
                      align_trim=(0.05,0.95),
                      hitlist_size=5000,
                      e_value_cutoff=0.001,
                      gapcosts=(11,1),
                      num_ncbi_threads=1,
                      num_recip_blast_threads=-1,
                      overwrite=False,
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
    target_seq_number : int, default=500
        desired number of sequences in final alignment
    ncbi_blast_db : str or None, default="nr"
        NCBI blast database to use. If None, use a local database. Incompatible
        with local_blast_db.
    local_blast_db : str, optional
        Local blast database to use. If None, use an NCBI database. Incompatible
        with ncbi_blast_db.
    local_recip_blast_db : str, optional
        Local blast database to use for reciprocal blast. If None, construct a
        reciprocal blast database by downloading the proteomes of the key
        species from the ncbi.
    within_species_redundancy_cutoff : float, default=0.99
        merge sequences with sequence identity above cutoff when removing
        redundancy of sequences within species.
    sparse_column_cutoff : float, default=0.95
        when checking alignment quality, a column is sparse if it has gaps in
        more than sparse_column_cutoff sequences.
    align_trim : tuple, default=(0.05,0.95)
        when checking alignment quality, do not score the first and last parts
        of the alignment. Interpreted like a slice, but with percentages.
        (0.0,1.0) would not trim; (0.05,0,98) would trim the first 0.05 off the
        front and the last 0.02 off the back.
    hitlist_size : int, default=5000
        download only the top hitlist_size hits
    e_value_cutoff : float, default=0.001
        only take hits with e_value better than e_value_cutoff
    gapcost : tuple, default=(11,1)
        BLAST gapcosts (length 2 tuple of ints)
    num_ncbi_threads : int, default=1
        number of threads to use for the NCBI blast. -1 means use all available.
        (Multithreading rarely speeds up remote BLAST).
    num_recip_blast_threads : int, default=-1
        number of threads to use for local reciprocal blast. -1 means all
        available.
    overwrite : bool, default=False
        overwrite out_dir if it already exists
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
    #
    # # Try to create the output directory
    # out_dir = str(out_dir)
    # if os.path.exists(out_dir):
    #     if os.path.isdir(out_dir):
    #         if overwrite:
    #             shutil.rmtree(out_dir)
    #         else:
    #             err = f"\nout_dir '{out_dir}' already exists\n\n"
    #             raise FileExistsError(err)
    #     else:
    #         err = f"\nout_dir '{out_dir}' already exists and is not a directory\n\n"
    #         raise FileExistsError(err)
    # os.mkdir(out_dir)
    #
    # If the seed_df is a file, copy that into the output directory
    if type(seed_df) is str:
        if os.path.exists(seed_df):
            root = f"00_{os.path.split(seed_df)[-1]}"
            shutil.copy(seed_df,os.path.join(out_dir,root))
            seed_df = root
        else:
            err = f"\nseed_df '{seed_df}' not found.\n\n"
            raise FileNotFoundError(err)

    # Change into the output directory
    cwd = os.getcwd()
    os.chdir(out_dir)

    step_counter = 1

    print("-------------------------------------------------------------------")
    print("Building initial topiary dataframe.")
    print("-------------------------------------------------------------------")
    print("",flush=True)

    kwargs = {"seed_df":seed_df,
              "ncbi_blast_db":ncbi_blast_db,
              "local_blast_db":local_blast_db,
              "hitlist_size":hitlist_size,
              "e_value_cutoff":e_value_cutoff,
              "gapcosts":gapcosts,
              "num_threads":num_ncbi_threads}

    out = topiary.df_from_seed(**kwargs)

    df = out[0]
    phylo_context = out[1]
    key_species = out[2]
    paralog_patterns = out[3]

    topiary.write_dataframe(df,f"{step_counter:02d}_initial-dataframe.csv")
    step_counter += 1

    print("-------------------------------------------------------------------")
    print("Doing reciprocal blast.")
    print("-------------------------------------------------------------------")
    print("",flush=True)

    # If no reciprocal blast database is specified, construct from key species. 
    if local_recip_blast_db is None:

        proteome_list = []
        for k in key_species:
           proteome_list.append(topiary.ncbi.get_proteome(species=k))

        base = "".join([random.choice(string.ascii_letters) for _ in range(10)])
        blast_db = f"{base}_local_blast"
        topiary.ncbi.make_blast_db(proteome_list,"blast_db",overwrite=True)

        local_recip_blast_db = "blast_db"

    df = topiary.recip_blast(df,
                             paralog_patterns=paralog_patterns,
                             local_blast_db="blast_db",
                             num_threads=num_recip_blast_threads)

    topiary.write_dataframe(df,f"{step_counter:02d}_recip-blast-dataframe.csv")
    step_counter += 1

    print("-------------------------------------------------------------------")
    print("Reducing number of sequences.")
    print("-------------------------------------------------------------------")
    print("",flush=True)

    kwargs = {"df":df,
              "paralog_column":"recip_paralog",
              "target_seq_number":target_seq_number,
              "key_species":key_species,
              "within_species_redundancy_cutoff":within_species_redundancy_cutoff,
              "sparse_column_cutoff":sparse_column_cutoff,
              "align_trim":align_trim,
              "verbose":verbose}

    df = topiary.taxonomic_sample(**kwargs)
    topiary.write_dataframe(df,f"{step_counter:02d}_sampled-dataframe.csv")
    step_counter += 1

    print("-------------------------------------------------------------------")
    print("Aligning sequences.")
    print("-------------------------------------------------------------------")
    print("",flush=True)

    df = topiary.run_muscle(df)
    topiary.write_dataframe(df,f"{step_counter:02d}_aligned-dataframe.csv")
    step_counter += 1

    topiary.write_fasta(df,f"{step_counter:02d}_alignment.fasta",seq_column="alignment")

    os.chdir(cwd)

    return df
