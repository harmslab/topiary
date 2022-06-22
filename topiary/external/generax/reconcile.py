"""
Reconcile a gene tree with a species tree using generax.
"""

import topiary
from topiary.external._interface import prep_calc, write_run_information

from ._generax import setup_generax, run_generax, GENERAX_BINARY

import ete3
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

    Parameters
    ----------
    previous_dir : str, optional
        directory containing previous calculation. function will grab the the
        csv, model, and tree from the previous run. If this is not specified,
        `df`, `model`, and `tree_file` arguments must be specified.
    df : pandas.DataFrame or str, optional
        topiary data frame or csv written out from topiary df. Will override
        dataframe from `previous_dir` if specified.
    model : str, optional
        model (i.e. "LG+G8"). Will override model from `previous_dir`
        if specified.
    tree_file : str
        tree_file in newick format. Will override tree from `previous_dir` if
        specified.
    allow_horizontal_transfer : bool, default=False
        whether to allow horizontal transfer during reconcilation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    output: str, optional
        output directory. If not specified, create an output directory with
        form "generax_reconcilation_randomletters"
    overwrite : bool, default=False
        whether or not to overwrite existing output directory
    generax_binary : str, optional
        what generax binary to use

    Returns
    -------
    Python.core.display.Image or None
        if running in jupyter notebook, return Image showing reconciled tree;
        otherwise, return None.
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


    # Get outgroups (e.g. leaves descending from each half after root)
    reconcile_file = os.path.join("working","result","reconciliations","reconcile_events.newick")
    reconcile_tree = ete3.Tree(reconcile_file,format=1)
    root = reconcile_tree.get_tree_root()
    root_children = root.get_children()
    outgroup = [[n.name for n in r.get_leaves()] for r in root_children]

    # Write run information
    write_run_information(outdir=outdir,
                          df=df,
                          calc_type="reconciliation",
                          model=model,
                          cmd=cmd,
                          outgroup=outgroup)

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
