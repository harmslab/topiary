"""
Generate maximum likelihood tree from an alignment given an evolutionary
model.
"""

import topiary

from ._raxml import run_raxml
from ._raxml import RAXML_BINARY
from topiary._private.supervisor import Supervisor

import os, shutil, glob

def generate_ml_tree(previous_dir=None,
                     df=None,
                     model=None,
                     gene_tree=None,
                     calc_dir="ml_tree",
                     overwrite=False,
                     bootstrap=False,
                     supervisor=None,
                     num_threads=-1,
                     raxml_binary=RAXML_BINARY):
    """
    Generate maximum likelihood tree from an alignment given an evolutionary
    model.

    Parameters
    ----------
    previous_dir : str, optional
        directory containing previous calculation. function will grab the the
        csv, model, and tree from the previous run. If this is not specified,
        `df`, `model`, and `gene_tree` arguments must be specified.
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
    calc_dir : str, default="ml_tree"
        calculation directory. Will be created.
    overwrite : bool, default=False
        whether or not to overwrite existing output
    bootstrap : bool, default=False
        whether or not to do bootstrap replicates
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

    calc_type = "ml_tree"
    if bootstrap:
        calc_type = "ml_bootstrap"

    supervisor.create_calc_dir(calc_dir=calc_dir,
                               calc_type=calc_type,
                               overwrite=overwrite,
                               df=df,
                               gene_tree=gene_tree,
                               model=model)

    supervisor.check_required(required_values=["model"],
                              required_files=["alignment.phy","dataframe.csv"])

    os.chdir(supervisor.working_dir)

    other_args = []

    # If we're doing bootstrapping
    if bootstrap:
        algorithm = "--all"
        other_args.extend(["--bs-trees","autoMRE","--bs-write-msa"])
    else:
        algorithm = "--search"

    # Run raxml to create tree
    cmd = run_raxml(run_directory="infer-ml-tree",
                    algorithm=algorithm,
                    alignment_file=supervisor.alignment,
                    tree_file=supervisor.gene_tree,
                    model=supervisor.model,
                    seed=supervisor.seed,
                    supervisor=supervisor,
                    num_threads=num_threads,
                    raxml_binary=raxml_binary,
                    other_args=other_args)


    # Grab the final tree and store as tree.newick
    supervisor.stash(os.path.join("infer-ml-tree","alignment.phy.raxml.bestTree"),
                     target_name="gene-tree.newick")

    # If we ran bootstrap, get tree with supports and bootstrap information
    if bootstrap:

        supervisor.stash(os.path.join("infer-ml-tree","alignment.phy.raxml.support"),
                         target_name="gene-tree_supports.newick")

        bsmsa = glob.glob(os.path.join("infer-ml-tree",
                                       "alignment.phy.raxml.bootstrapMSA.*.phy"))
        for b in bsmsa:
            number = int(b.split(".")[-2])
            target_name = os.path.join("bootstrap_replicates",
                                       f"bsmsa_{number:04d}.phy")
            supervisor.stash(b,target_name=target_name)

        supervisor.stash(os.path.join("infer-ml-tree",
                                      "alignment.phy.raxml.bootstraps"),
                         target_name=os.path.join("bootstrap_replicates",
                                                  "bootstraps.newick"))

    # Close out
    os.chdir(supervisor.starting_dir)
    return supervisor.finalize(successful=True,plot_if_success=True)
