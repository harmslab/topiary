__description__ = \
"""
Reconcile a gene tree with a species tree using generax.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-05-16"

import topiary
from topiary.external.interface import prep_calc, write_run_information

from ._generax import setup_generax, run_generax, GENERAX_BINARY

import numpy as np

import os, glob, shutil

def reconcile(previous_dir=None,
              df=None,
              model=None,
              tree_file=None,
              allow_horizontal_transfer=False,
              output=None,
              overwrite=False,
              generax_binary=GENERAX_BINARY):
    """
    Reoncile the gene tree to the species tree using generax.

    previous_dir: directory containing previous calculation. prep_calc will
                  grab the the csv, model, and tree from the previous run.
    df: topiary data frame or csv written out from topiary df
    model: model (e.g. LG+G8).
    tree_file: tree_file in newick format.
    allow_horizontal_transfer: whether to allow horizontal transfer during
                               reconcilation. If True, use the UndatedDTL model.
                               If False, use the UndatedDL model.
    output: output directory. If not specified, create an output directory with
            form "generax_reconcilation_randomletters"
    overwrite: whether or not to overwrite existing output (default False)
    threads: number of threads to use XXX
    generax_binary: what generax binary to use
    """

    # Prepare for the calculation, loading in previous calculation and
    # combining with arguments as passed in.
    result = prep_calc(previous_dir=previous_dir,
                       df=df,
                       model=model,
                       tree_file=model,
                       output=output,
                       overwrite=overwrite,
                       output_base="generax_reconcilation")

    df = result["df"]
    csv_file = result["csv_file"]
    model = result["model"]
    tree_file = result["tree_file"]
    alignment_file = result["alignment_file"]
    starting_dir = result["starting_dir"]

    required = [df,model,tree_file]
    for r in required:
        if r is None:
            err = "\nA dataframe, model, and tree are required for this "
            err += "calculation.\n\n"
            raise ValueError(err)

    # Set up generax directory
    setup_generax(df,tree_file,model,"working")

    # Actually run generax
    cmd = run_generax(run_directory="working",
                      allow_horizontal_transfer=allow_horizontal_transfer,
                      generax_binary=generax_binary)

    # Make output directory to hold final outputs
    outdir = "output"
    os.mkdir(outdir)

    # Copy in tree.newick
    shutil.copy(os.path.join("working","result","results","reconcile","geneTree.newick"),
                os.path.join("output","tree.newick"))

    # Write run information
    write_run_information(outdir=outdir,
                          df=df,
                          calc_type="reconciliation",
                          model=model,
                          cmd=cmd)

    # Copy reconcilation information
    shutil.copytree(os.path.join("working","result","reconciliations"),
                    os.path.join("output","reconcilations"))

    print(f"\nWrote results to {os.path.abspath(outdir)}\n")

    # Leave working directory
    os.chdir(starting_dir)

    # Write out a summary tree.
    ret = topiary.draw.reconciliation_tree(run_dir=output,
                                           output_file=os.path.join(output,
                                                                    "output",
                                                                    "summary-tree.pdf"))
    if topiary._in_notebook:
        return ret
