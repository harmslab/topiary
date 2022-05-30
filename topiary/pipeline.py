
import topiary
from topiary import _arg_processors
import pandas as pd
import numpy as np

import os, re

def read_paralog_patterns(paralog_patterns):

    err_summary = ["This file should have the following format:\n",
                   "paralog1;pattern1;pattern2;pattern3"
                   "paralog2;pattern4;pattern5;pattern6",
                   "...\n",
                   "For a real pair of proteins, this could be:\n",
                   "LY96;MD2;MD-2;ESOP1",
                   "LY86;MD1;MD-1\n",
                   "Blank lines and anything after # on a line are ignored\n"]
    err_summary = "\n".join(err_summary)

    if type(paralog_patterns) is str:
        if os.path.isfile(paralog_patterns):

            out_dict = {}
            with open(paralog_patterns) as f:
                for line in f:

                    # Remove everything after '#' (comment)
                    l = line.split("#")[0].strip()

                    # Skip blank lines
                    if l.strip() == "":
                        continue

                    # Now split line on ";" and get non-"" entries.
                    col = [c.strip() for c in l.split(";")]
                    col = [c for c in col if c != ""]

                    if len(col) < 2:
                        err = f"\nCould not parse line '{line}' in paralog_patterns.\n\n"
                        err += err_summary
                        raise ValueError(err)

                    # Create pattern dictionary, searching for all unique patterns
                    out_dict[col[0]] = list(set(col))

            paralog_patterns = out_dict

    patterns = _arg_processors.process_paralog_patterns(paralog_patterns)

    return dict([(p[1],p[0]) for p in patterns])



def rockit(input,
           paralog_patterns,
           output_dir="rockit",
           redundancy_cutoff=0.85,
           target_number_of_sequences=500,
           ncbi_rev_blast_db="nr",
           local_rev_blast_db=None,
           ncbi_rev_blast_taxid=9606,
           phylo_context="All life"):
    """
    """
    step_counter = 0

    # Read the dictionary or file containing paralog_patterns
    paralog_patterns = read_paralog_patterns(paralog_patterns)

    # Read the dataframe containing blast xml files
    df = topiary.ncbi_blast_xml_to_df(input)
    topiary.write_dataframe(df,"{step_counter:02d}_initial-dataframe.csv")
    step_counter += 1

    # Assign nicknames to all sequences using paralog patterns
    df = topiary.create_nicknames(df,paralog_patterns=paralog_patterns)
    topiary.write_dataframe(df,"{step_counter:02d}_nickname-dataframe.csv")
    step_counter += 1

    # Get ott id. We do this early on in case we can't find a species ott.
    df = topiary.get_ott_id(df,phylo_context=phylo_context)
    topiary.write_dataframe(df,"{step_counter:02d}_ott-dataframe.csv")
    step_counter += 1

    # Now remove redundancy.
    df = topiary.remove_redundancy(df,cutoff=redundancy_cutoff)

    # If a reverse blast database got passed in, do reerse blast.
    if local_rev_blast_db is not None or ncbi_rev_blast_db is not None:

        df = topiary.reverse_blast(df,
                                   paralog_patterns=paralog_patterns,
                                   local_rev_blast_db=local_rev_blast_db,
                                   ncbi_rev_blast_db=ncbi_rev_blast_db,
                                   taxid=ncbi_rev_blast_taxid)

    # Align sequences
    df = topiary.run_muscle(df)

    # Clean alignment
    df = topiary.quality.clean_alignment(df)

    # Align final set of high quality sequences
    df = topiary.run_muscle(df)

    topiary.write_fasta(df,"big-alignment.fasta",seq_column="alignment",overwrite=True)

def generate_ancestors(df,output_dir,do_bootstrap=True):
    """
    """

    topiary.find_best_model(df,
                            output="00_find-model")

    topiary.generate_ml_tree(previous_dir="00_find-model",
                             output="01_ml-tree",
                             bootstrap=do_bootstrap)

    topiary.reconcile(previous_dir="01_ml-tree",
                      output="02_reconciliation")

    topiary.generate_ancestors(previous_dir="02_reconciliation",
                               output="03_ancestors")
