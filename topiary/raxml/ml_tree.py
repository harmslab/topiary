__description__ = \
"""
Generate maximum likelihood ancestors using raxml.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2021-07-22"

import topiary

from ._raxml import prep_calc, run_raxml, RAXML_BINARY

import os, shutil

def generate_ml_tree(df,
                     model,
                     tree_file=None,
                     output=None,
                     threads=1,
                     raxml_binary=RAXML_BINARY,
                     bootstrap=True,
                     write_bs_msa=True):
    """
    Generate maximum likelihood tree from an alignment given a substitution
    model.

    df: topiary data frame or csv written out from topiary df
    model: model (e.g. LG+G8).
    tree_file: tree_file in newick format. If not specified, a parsimony tree
               will be generated. used as starting point.
    output: name of output directory.
    threads: number of threads to use
    raxml_binary: what raxml binary to use
    bootstrap: whether or not to do bootstrap replicates
    write_bs_msa: whether or not to write out all bootstrap alignments
    """

    # Copy files in, write out alignment, move into working directory, etc.
    result = prep_calc(df=df,
                       output=output,
                       other_files=[tree_file],
                       output_base="generate_ml_tree_")

    df = result["df"]
    csv_file = result["csv_file"]
    alignment_file = result["alignment_file"]
    tree_file = result["other_files"][0]
    starting_dir = result["starting_dir"]

    other_args = []

    # If we're doing bootstrapping
    if bootstrap:

        # Configure bootstrap
        other_args.extend(["--bs-trees","autoMRE"])

        # write out bootstrap alignments
        if write_bs_msa:
            other_args.append("--bs-write-msa")

    # Run raxml to create tree
    run_raxml(algorithm="--all",
              alignment_file=alignment_file,
              tree_file=tree_file,
              model=model,
              dir_name="01_make-ml-tree",
              seed=True,
              threads=threads,
              raxml_binary=raxml_binary,
              other_args=other_args)

    # Write out a pretty version of the tree
    shutil.copy("01_make-ml-tree/alignment.raxml.support","02_ml-tree.newick")

    # Leave working directory
    os.chdir(starting_dir)
