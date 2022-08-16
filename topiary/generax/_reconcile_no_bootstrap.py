"""
Reconcile a gene tree with a species tree using generax without bootstraps.
"""

import topiary

from topiary._private import check
from ._generax import setup_generax
from ._generax import run_generax
from ._generax import GENERAX_BINARY

import ete3
import numpy as np

import os, glob, shutil


def reconcile_no_bootstrap(df,
                           model,
                           gene_tree,
                           species_tree,
                           allow_horizontal_transfer,
                           seed,
                           overwrite,
                           supervisor,
                           num_threads,
                           generax_binary):
    """
    Reconcile the gene tree to the species tree using generax. This should be
    called by generax.reconcile.reconcile
    Parameters
    ----------
    df : pandas.DataFrame or str, optional
        topiary data frame or csv written out from topiary df. Will override
        dataframe from `prev_calculation` if specified.
    model : str, optional
        model (i.e. "LG+G8"). Will override model from `prev_calculation`
        if specified.
    gene_tree : str, optional
        gene_tree in newick format.
    species_tree : str, optional
        species tree in newick format.
    allow_horizontal_transfer : bool, default=True
        whether to allow horizontal transfer during reconcilation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    seed : bool,int,str
        If true, pass a randomly generated seed to raxml. If int or str, use
        that as the seed. (passed via --seed)
    overwrite : bool
        whether or not to overwrite existing output directory
    supervisor : Supervisor
        supervisor instance to keep track of inputs and outputs
    num_threads : int
        number of threads to use. if -1 use all available.
    generax_binary : str
        what generax binary to use

    Returns
    -------
    plot : toyplot.canvas or None
        if running in jupyter notebook, return toyplot.canvas; otherwise, return
        None.
    """

    os.chdir(supervisor.working_dir)

    supervisor.event("Setting up reconcilation directory")

    # Set up generax directory
    setup_generax(df=df,
                  gene_tree=gene_tree,
                  model=model,
                  out_dir="generax",
                  species_tree=species_tree)

    # Actually run generax
    cmd = run_generax(run_directory="generax",
                      allow_horizontal_transfer=allow_horizontal_transfer,
                      seed=seed,
                      supervisor=supervisor,
                      num_threads=num_threads,
                      generax_binary=generax_binary)

    # Get species tree used for this calculation
    supervisor.stash(os.path.join(supervisor.working_dir,"generax","species_tree.newick"),"species-tree.newick")

    # Copy in tree.newick
    supervisor.stash(os.path.join(supervisor.working_dir,
                                  "generax",
                                  "result",
                                  "results",
                                  "reconcile",
                                  "geneTree.newick"),
                     "reconciled-tree.newick")

    # Copy reconcilation information
    reconcile_dir = os.path.join(supervisor.working_dir,
                                 "generax",
                                 "result",
                                 "reconciliations")
    supervisor.stash(os.path.join(reconcile_dir,"reconcile_events.newick"),
                     "reconciled-tree_events.newick")
    supervisor.stash(reconcile_dir,"reconcilations")

    os.chdir(supervisor.starting_dir)
    return supervisor.finalize(successful=True,plot_if_success=True)
