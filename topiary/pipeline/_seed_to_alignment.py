__description__ = \
"""
Pipeline that takes a seed dataframe, BLASTS to find sequence hits,
performs quality control, lowers alignment redundancy in a taxonomically
informed fashion, and generates an alignment.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-17"

import topiary

import numpy as np
import pandas as pd

import random, string


def seed_to_alignment(seed_df,
                      target_seq_number=500,
                      ncbi_blast_db="nr",
                      local_blast_db=None,
                      within_species_redundancy_cutoff=0.99,
                      sparse_column_cutoff=0.95,
                      sparse_run_length_keep_percentile=0.98,
                      fx_missing_dense_cutoff=0.9,
                      align_trim=(0.05,0.95),
                      hitlist_size=5000,
                      e_value_cutoff=0.001,
                      gapcosts=(11,1),
                      num_ncbi_threads=1,
                      num_recip_blast_threads=-1,
                      verbose=False):
    """
    Pipeline that takes a seed dataframe, BLASTS to find sequence hits,
    performs quality control, lowers alignment redundancy in a taxonomically
    informed fashion, and generates an alignment.

    Parameters
    ----------
        seed_df: dataframe with at least four columns: name, species, sequence,
            and aliases. See documentation on seed dataframes for details.
        target_seq_number: desired number of sequences in final alignment
        ncbi_blast_db: database on ncbi against which to blast (incompatible
            with local_blast_db)
        local_blast_db: local database against which to blast (incompatible
            with ncbi_blast_db)
        within_species_redundancy_cutoff: merge sequences with sequence identity
            above cutoff when removing redundancy of sequences within species.
        sparse_column_cutoff: when checking alignment quality, a column is
            sparse if it has gaps in more than sparse_column_cutoff sequences.
        sparse_run_length_keep_percentile: when checking alignment quality,
            remove the sequences with the longest insertions. Toss the longest
            (1 - sparse_run_length_keep_percentile)*num_sequences sequences.
        fx_missing_dense_cutoff: when checking alignment quality, remove
            sequences that are missing more than 1 - fx_missing_dense_cutoff
            of the columns that are not sparse.
        align_trim: when checking alignment quality, do not score the first and
            last parts of the alignment. Interpreted like a slice, but with
            percentages. (0.0,1.0) would not trim; (0.05,0,98) would trim the
            first 0.05 off the front and the last 0.02 off the back.
        hitlist_size: download only the top hitlist_size hits
        e_value_cutoff: only take hits with e_value better than e_value_cutoff
        gapcost: BLAST gapcosts (length 2 tuple of ints)
        num_ncbi_threads: number of threads to use for the NCBI blast. -1 means
            use all available. (Multithreading rarely speeds up remote BLAST).
        num_recip_blast_threads: number of threads to use for local reciprocal
            blast. -1 means all available.
        verbose: verbosity of output

    Return
    ------
        Topiary dataframe with aligned, quality-controlled sequences.
    """

    step_counter = 1

    print("-------------------------------------------------------------------")
    print("Building initial topiary dataframe.")
    print("-------------------------------------------------------------------")
    print("",flush=True)

    kwargs = {"df":seed_df,
              "ncbi_blast_db":ncbi_blast_db,
              "local_blast_db":local_blast_db,
              "hitlist_size":hitlist_size,
              "e_value_cutoff":e_value_cutoff,
              "gapcost":gapcosts,
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

    proteome_list = []
    for k in key_species:
        proteome_list.append(topiary.ncbi.get_proteome(species=k))

    base = "".join([random.choice(string.ascii_letters) for _ in range(10)])
    blast_db = f"{base}_local_blast"
    topiary.ncbi.make_blast_db(proteome_list,"blast_db",overwrite=True)
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
              "sparse_run_length_keep_percentile":sparse_run_length_keep_percentile,
              "fx_missing_dense_cutoff":fx_missing_dense_cutoff,
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

    return df
