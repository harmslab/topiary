"""
Reconcile a gene tree with a species tree using generax.
"""

import topiary
from topiary._private.interface import prep_calc, write_run_information
from topiary._private.threads import get_num_threads

from topiary._private import check
from ._generax import setup_generax, run_generax, GENERAX_BINARY
from ._reconcile_bootstrap import reconcile_bootstrap

import ete3
import numpy as np

import os, glob, shutil

def reconcile(previous_dir=None,
              df=None,
              model=None,
              tree_file=None,
              allow_horizontal_transfer=True,
              output=None,
              overwrite=False,
              bootstrap=False,
              use_mpi=True,
              num_threads=-1,
              num_cores=None,
              generax_binary=GENERAX_BINARY):
    """
    Reconcile the gene tree to the species tree using generax.

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
    allow_horizontal_transfer : bool, default=True
        whether to allow horizontal transfer during reconcilation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    output: str, optional
        output directory. If not specified, create an output directory with
        form "generax_reconcilation_randomletters"
    overwrite : bool, default=False
        whether or not to overwrite existing output directory
    bootstrap : bool, default=False
        whether or not to do bootstrap replicates. Requires a previous_dir from
        an ML tree estimation that was run with bootstrap replicates.
    use_mpi : bool, default=True
        whether or not to use mpi
    num_threads : int, default=-1
        number of threads to use. if -1 use all available.
    num_cores : int, optional
        number of cores. if None, infer. (Useful for mpirun calls).
    generax_binary : str, optional
        what generax binary to use

    Returns
    -------
    plot : toyplot.canvas or None
        if running in jupyter notebook, return toyplot.canvas; otherwise, return
        None.
    """

    # If bootstrap requested, run bootstrap version
    bootstrap = check.check_bool(bootstrap,"bootstrap")
    if bootstrap:
        ret = reconcile_bootstrap(previous_dir=previous_dir,
                                  df=df,
                                  model=model,
                                  tree_file=tree_file,
                                  allow_horizontal_transfer=allow_horizontal_transfer,
                                  output=output,
                                  overwrite=overwrite,
                                  use_mpi=use_mpi,
                                  num_threads=num_threads,
                                  num_cores=num_cores,
                                  generax_binary=GENERAX_BINARY)
        return ret


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
    output = result["output"]
    existing_trees = result["existing_trees"]

    required = [df,model,tree_file]
    for r in required:
        if r is None:
            err = "\nA dataframe, model, and tree are required for this "
            err += "calculation.\n\n"
            raise ValueError(err)

    # Set up generax directory
    setup_generax(df,tree_file,model,"working")

    # Get number of threads to use. If no mpi, use only one thread
    num_threads = get_num_threads(num_threads,num_cores)
    if not use_mpi:
        num_threads = 1

    # Actually run generax
    cmd = run_generax(run_directory="working",
                      allow_horizontal_transfer=allow_horizontal_transfer,
                      num_threads=num_threads,
                      generax_binary=generax_binary)

    # Make output directory to hold final outputs
    os.mkdir("output")

    # Copy trees from previous calculation in. This will preserve any that our
    # new calculation did not wipe out.
    for t in existing_trees:
        tree_filename = os.path.split(t)[-1]
        shutil.copy(t,os.path.join("output",tree_filename))

    # Copy in tree.newick
    shutil.copy(os.path.join("working","result","results","reconcile","geneTree.newick"),
                os.path.join("output","tree.newick"))

    # Copy reconcilation information
    shutil.copytree(os.path.join("working","result","reconciliations"),
                    os.path.join("output","reconcilations"))
    shutil.copy(os.path.join("output","reconcilations","reconcile_events.newick"),
                os.path.join("output","tree_events.newick"))

    # Write run information
    write_run_information(outdir="output",
                          df=df,
                          calc_type="reconciliation",
                          model=model,
                          cmd=cmd)

    print(f"\nWrote results to {os.path.abspath('output')}\n")

    # Leave working directory
    os.chdir(starting_dir)

    # Write out a summary tree.
    return topiary.draw.tree(run_dir=output,
                             output_file=os.path.join(output,
                                                      "output",
                                                      "summary-tree.pdf"))
