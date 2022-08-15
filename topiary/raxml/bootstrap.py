"""
Generate bootstrap replicates for an existing tree and then calculate boostrap
supports.
"""

import topiary

from ._raxml import run_raxml, RAXML_BINARY
from topiary._private.supervisor import Supervisor

import os
import shutil
import glob

def generate_bootstraps(previous_dir=None,
                        df=None,
                        model=None,
                        gene_tree=None,
                        calc_dir="ml_bootstrap",
                        overwrite=False,
                        num_bootstraps=None,
                        supervisor=None,
                        num_threads=-1,
                        raxml_binary=RAXML_BINARY):
    """
    Generate bootstrap replicates for an existing tree and then calculate boostrap
    supports.

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
    gene_tree : str or ete3.Tree or dendropy.Tree
        gene_tree. Used as starting point for calculation. Will override tree
        from `previous_dir` if specified. Should be newick with only leaf names
        and branch lengths. If this an ete3 or dendropy tree, it will be written
        out with leaf names and branch lengths; all other data will be dropped
    calc_dir : str, default="ml_bootstrap"
        directory in which to do calculation.
    overwrite : bool, default=False
        whether or not to overwrite existing calc_dir
    num_bootstraps : int, optional
        how many bootstrap replicates to generate. If None, use autoMRE to
        automatically infer the number of replicates given the data.
    supervisor : Supervisor, optional
        instance of Supervisor for managing calculation inputs and outputs
    num_threads : int, default=-1
        number of threads to use. if -1, use all avaialable
    raxml_binary : str, optional
        what raxml binary to use

    Returns
    -------
    plot : toyplot.canvas or None
        if running in jupyter notebook, return toyplot.canvas; otherwise, return
        None.
    """

    if supervisor is None:
        supervisor = Supervisor(calc_dir=previous_dir)

    supervisor.create_calc_dir(calc_dir=calc_dir,
                               calc_type="ml_bootstrap",
                               overwrite=overwrite,
                               df=df,
                               gene_tree=gene_tree,
                               model=model)

    supervisor.check_required(required_values=["model"],
                              required_files=["alignment.phy","dataframe.csv",
                                              "gene-tree.newick"])

    os.chdir(supervisor.working_dir)

    # Figure out how to do bootstrapping
    algorithm = "--bootstrap"
    if num_bootstraps is None:
        num_bootstraps = "autoMRE"
    else:
        num_bootstraps = f"{num_bootstraps:d}"

    other_args = []
    other_args.extend(["--bs-trees",num_bootstraps,"--bs-write-msa"])

    # Run raxml to generate bootstrap replicates
    cmd1 = run_raxml(run_directory="00_generate-replicates",
                     algorithm=algorithm,
                     alignment_file=supervisor.alignment,
                     tree_file=supervisor.gene_tree,
                     model=supervisor.model,
                     seed=supervisor.seed,
                     other_args=other_args,
                     supervisor=supervisor,
                     num_threads=num_threads,
                     raxml_binary=raxml_binary)

    # Get bs_file
    bs_file = os.path.abspath(os.path.join("00_generate-replicates",
                                           "alignment.phy.raxml.bootstraps"))

    cmd2 = run_raxml(run_directory="01_combine-bootstraps",
                     algorithm="--support",
                     tree_file=supervisor.gene_tree,
                     log_to_stdout=False,
                     suppress_output=True,
                     other_args=["--bs-trees","alignment.phy.raxml.bootstraps","--redo"],
                     other_files=[bs_file],
                     supervisor=supervisor,
                     num_threads=1,
                     raxml_binary=raxml_binary)


    # Get newick files from previous output directory and put in new output
    supervisor.copy_output_to_output("*.newick")
    supervisor.stash(os.path.join(supervisor.input_dir,"dataframe.csv"))

    # Grab tree with bootstraps and store as tree_supports.newick
    supervisor.stash(os.path.join("01_combine-bootstraps",
                                  "tree.newick.raxml.support"),
                     "gene-tree_supports.newick")

    # Copy bootstrap results to the output directory
    bs_out = "bootstrap_replicates"
    bsmsa = glob.glob(os.path.join("00_generate-replicates",
                                   "alignment.phy.raxml.bootstrapMSA.*.phy"))
    for b in bsmsa:
        number = int(b.split(".")[-2])
        supervisor.stash(b,os.path.join(bs_out,f"bsmsa_{number:04d}.phy"))

    supervisor.stash(os.path.join("01_combine-bootstraps",
                                  "alignment.phy.raxml.bootstraps"),
                     os.path.join(bs_out,"bs-trees.newick"))

    os.chdir(supervisor.starting_dir)
    return supervisor.finalize(successful=True,plot_if_success=True)
