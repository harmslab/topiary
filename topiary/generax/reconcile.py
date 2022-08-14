"""
Reconcile a gene tree with a species tree using generax.
"""

import topiary
from topiary._private import Supervisor
from topiary._private.mpi import get_num_slots
from topiary._private.mpi import check_mpi_configuration

from ._reconcile_bootstrap import reconcile_bootstrap
from ._reconcile_no_bootstrap import reconcile_no_bootstrap

from ._generax import GENERAX_BINARY

import subprocess
import os

def reconcile(previous_dir=None,
              df=None,
              model=None,
              tree_file=None,
              species_tree_file=None,
              allow_horizontal_transfer=True,
              bootstrap=False,
              calc_dir="generax_reconcilation",
              overwrite=False,
              supervisor=None,
              num_threads=-1,
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
    tree_file : str, optional
        tree_file in newick format. Will override tree from `previous_dir` if
        specified.
    species_tree_file : str, optional
        species tree in newick format.
    allow_horizontal_transfer : bool, default=True
        whether to allow horizontal transfer during reconcilation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    bootstrap: bool, default=False
        whether or not to do bootstrap replicates. if True, previous_dir must
        point to a raxml ml_bootstrap run
    calc_dir: str, default="generax_reconcilation"
        name of calc_dir directory
    overwrite : bool, default=False
        whether or not to overwrite existing calc_dir directory
    supervisor : Supervisor, optional
        supervisor instance to keep track of inputs and outputs
    num_threads : int, default=-1
        number of threads to use. if -1 use all available.
    generax_binary : str, optional
        what generax binary to use

    Returns
    -------
    plot : toyplot.canvas or None
        if running in jupyter notebook, return toyplot.canvas; otherwise, return
        None.
    """

    # Make sure that generax is in the path
    try:
        subprocess.run([generax_binary],capture_output=True)
    except FileNotFoundError:
        err = f"\ngenerax binary '{generax_binary}' not found in path\n\n"
        raise ValueError(err)

    # Get number of slots
    if num_threads == -1:
        num_threads = get_num_slots(generax_binary)

    # Check sanity of mpi configuration/number of slots
    check_mpi_configuration(num_threads,generax_binary)

    # --------------------------------------------------------------------------
    # Load/parse calculation inputs

    # Load existing
    if supervisor is None:
        supervisor = Supervisor(previous_dir)

    # Create a calculation directory
    supervisor.create_calc_dir(calc_dir=calc_dir,
                               calc_type="reconcile",
                               overwrite=overwrite,
                               df=df,
                               tree=tree_file)

    if species_tree_file is not None:
        supervisor.stash(species_tree_file,
                         "species_tree.newick",
                         target_dir="input")

    if model is not None:
        supervisor.update("model",str(model))

    supervisor.update("allow_horizontal_transfer",allow_horizontal_transfer)

    supervisor.check_required(required_values=["model","allow_horizontal_transfer"],
                              required_files=["alignment.phy","dataframe.csv",
                                              "tree.newick"])

    if not bootstrap:

        return reconcile_no_bootstrap(previous_dir=previous_dir,
                                      df=df,
                                      model=model,
                                      tree_file=tree_file,
                                      species_tree_file=species_tree_file,
                                      allow_horizontal_transfer=allow_horizontal_transfer,
                                      calc_dir=calc_dir,
                                      overwrite=overwrite,
                                      supervisor=supervisor,
                                      num_threads=num_threads,
                                      generax_binary=generax_binary)

    else:

        try:
            prev_calc_type = supervisor.previous_entries[-1]["calc_type"]
        except:
            prev_calc_type = None

        if prev_calc_type != "ml_bootstrap":
            err = "\nboostrap reconciliation can only be started using a\n"
            err += "previous 'ml_bootstrap' run as its input.\n"
            raise ValueError(err)

        # Make sure bootstrap directory exists
        bs_dir = os.path.join(supervisor.previous_entries[-1]["output_dir"],
                              "bootstrap_replicates")

        if not os.path.isdir(bs_dir):
            err = f"\ninput directory '{supervisor.previous_entries[-1]['output_dir']}'\n"
            err += "does not have an bootstrap_replicates directory. Was\n"
            err += "this calculation run with bootstrap=True?\n\n"
            raise FileNotFoundError(err)

        return reconcile_bootstrap(supervisor=supervisor,
                                   calc_dir=calc_dir,
                                   overwrite=overwrite,
                                   num_threads=num_threads,
                                   generax_binary=generax_binary)
