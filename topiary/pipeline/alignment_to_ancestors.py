"""
Pipeline that starts from an alignment, finds the best phylogenetic model,
builds a maximum likelihood tree, reconciles this tree with the species tree,
and then infers ancestral proteins.
"""

import topiary
from topiary.raxml import RAXML_BINARY
from topiary.generax import GENERAX_BINARY
from topiary._private import installed
from topiary._private import software_requirements
from topiary._private import check
from topiary._private.mpi import check_mpi_configuration
from topiary.reports import create_report

import os
import random
import string
import shutil

def _check_restart(output,restart):

    run_calc = True
    if restart:

        # See if json file is there. If so, the run is done.
        json_file = os.path.join(output,"run_parameters.json")
        if os.path.isfile(json_file):
            run_calc = False
        else:
            # Nuke partial directory
            if os.path.isdir(output):
                shutil.rmtree(output)

    return run_calc

def alignment_to_ancestors(df,
                           out_dir=None,
                           starting_tree=None,
                           no_bootstrap=False,
                           force_reconcile=False,
                           force_no_reconcile=False,
                           horizontal_transfer=False,
                           alt_cutoff=0.25,
                           model_matrices=["cpREV","Dayhoff","DCMut","DEN","Blosum62",
                                           "FLU","HIVb","HIVw","JTT","JTT-DCMut","LG",
                                           "mtART","mtMAM","mtREV","mtZOA","PMB",
                                           "rtREV","stmtREV","VT","WAG"],
                           model_rates=["","G8"],
                           model_freqs=["","FC","FO"],
                           model_invariant=["","IC","IO"],
                           restart=False,
                           overwrite=False,
                           num_threads=-1,
                           raxml_binary=RAXML_BINARY,
                           generax_binary=GENERAX_BINARY):
    """
    Given an alignment, find the best phylogenetic model, build a maximum-
    likelihood tree, reconcile this tree with the species tree, and then infer
    ancestral protein sequences. Reconciliation is not done for df that only have
    bacterial sequences. User can force reconciliation to happen or not, 
    regardless of species in df, using force_reconcile and force_no_reconcile.

    Parameters
    ----------
    df : pandas.DataFrame or str
        topiary data frame or csv written out from topiary df.
    out_dir : str, optional
        output directory. If not specified, create an output directory with the
        format "alignment_to_ancestors_{randomletters}"
    starting_tree : str, optional
        tree in newick format. This will be used for the best model
        inference and starting tree for the maximum likelihood tree estimation.
        If not specified, the maximum parsimony tree is generated and used.
    no_bootstrap : bool, default=False
        do not do bootstrap replicates
    force_reconcile : bool, default=False
        reconcile gene and species trees
    force_no_reconcile : bool, default=False
        do not reconcile gene and species trees
    horizontal_transfer : bool, default=False
        whether to allow horizontal transfer during reconciliation. Default is
        to not allow transfer (UndatedDL). If set, use undated DTL. 
    alt_cutoff : float, default=0.25
        cutoff to use for altAll alternate ancestral protein sequence
        generation. Should be between 0 and 1.
    model_matrices : list, default=["cpREV","Dayhoff","DCMut","DEN","Blosum62","FLU","HIVb","HIVw","JTT","JTT-DCMut","LG","mtART","mtMAM","mtREV","mtZOA","PMB","rtREV","stmtREV","VT","WAG"]
        list of model matrices to check. If calling from command line, these
        can be specified directly (:code:`--model_matrices LG ...`) or by specifying
        a file with models on each line (:code:`--model_matrices SOME_FILE`)
    model_rates : list, default=["","G8"]
        ways to treat model rates. If calling from command line, these
        can be specified directly (:code:`--model_rates G8 ...`) or by specifying
        a file with rates on each line (:code:`--model_rates SOME_FILE`)
    model_freqs : list, default=["","FC","FO"]
        ways to treat model freqs. If calling from command line, these
        can be specified directly (:code:`--model_freqs FC FO ...`) or by specifying
        a file with freqs on each line (:code:`--model_freqs SOME_FILE`)
    model_invariant : list, default=["","IC","IO"]
        ways to treat invariant alignment columns. If calling from command line, these
        can be specified directly (:code:`--model_invariant IC IO ...`) or by specifying
        a file with invariants on each line (:code:`--model_invariant SOME_FILE`)
    restart : bool, default=False
        restart job from where it stopped in output directory. incompatible with
        overwrite
    overwrite : bool, default=False
        whether or not to overwrite existing output. incompatible with restart
    num_threads : int, default=-1
        number of threads to use. if -1 use all available
    raxml_binary : str, optional
        raxml binary to use
    generax_binary : str, optional
        what generax binary to use
    """

    # Read dataframe if string.
    if issubclass(type(df),str):
        df = topiary.read_dataframe(df)

    # Validate dataframe
    df = check.check_topiary_dataframe(df)
    
    # Validate starting_tree
    if starting_tree is not None:
        starting_tree = str(starting_tree)
        if not os.path.isfile(starting_tree):
            err = f"starting_tree '{starting_tree}' not found.\n"
            raise FileNotFoundError(err)

    # --------------------------------------------------------------------------
    # Flip logic from user interface (where flags turn off bootstrap and
    # reconciliation) to more readable flags that turn on (do_bootstrap,
    # do_reconcile)

    no_bootstrap = check.check_bool(no_bootstrap,"no_bootstrap")
    if no_bootstrap:
        do_bootstrap = False
    else:
        do_bootstrap = True

    # --------------------------------------------------------------------------
    # Decide how to do reconciliation based on user flags and taxonomic 
    # distribution. 

    force_no_reconcile = check.check_bool(force_no_reconcile,"force_no_reconcile")
    force_reconcile = check.check_bool(force_reconcile,"force_reconcile")
    if force_no_reconcile and force_reconcile:
        err = "force_no_reconcile and force_reconcile cannot both be set to True\n"
        raise ValueError(err)

    if force_reconcile:
        do_reconcile = True
    elif force_no_reconcile:
        do_reconcile = False
    else:

        # Figure out if the default is to reconcile or not based on taxonomic 
        # distribution of alignment. If only Bacteria or only Archaea, do not
        # try to reconcile. 
        mrca = topiary.opentree.ott_to_mrca(ott_list=list(df.ott),
                                            avoid_all_life=True)

        if mrca["is_microbial"]:
            do_reconcile = False
        else:
            do_reconcile = True

    # --------------------------------------------------------------------------
    # Validate calculation arguments

    # Convert to allow_horizontal_transfer
    horizontal_transfer = check.check_bool(horizontal_transfer,
                                           "horizontal_transfer")

    # alt-all cutoff
    alt_cutoff = check.check_float(alt_cutoff,
                                   "alt_cutoff",
                                   minimum_allowed=0,
                                   maximum_allowed=1)

    # model_matrices, model_freqs, model_rates, model_invariant go into
    # the first calculation (find_best_model) and have complicated validation.
    # Rely on that code to check.

    # --------------------------------------------------------------------------
    # Check sanity of overwrite, restart, and combination

    overwrite = check.check_bool(overwrite,"overwrite")
    restart = check.check_bool(restart,"restart")

    if overwrite and restart:
        err = "overwrite and restart flags are incompatible.\n"
        raise ValueError(err)

    num_threads = check.check_int(num_threads,
                                  "num_threads",
                                  minimum_allowed=-1)
    if num_threads == 0:
        err = "num_threads should be -1 (use all available) or a integer > 0\n"
        err += "indicating the number of threads to use.\n"
        raise ValueError(err)

    # --------------------------------------------------------------------------
    # Validate software stack required for this pipeline

    to_validate = [{"program":"raxml-ng",
                    "binary":raxml_binary,
                    "min_version":software_requirements["raxml-ng"],
                    "must_pass":True}]

    if do_reconcile:

        to_validate.append({"program":"generax",
                            "binary":generax_binary,
                            "min_version":software_requirements["generax"],
                            "must_pass":True})
        to_validate.append({"program":"mpirun",
                            "min_version":software_requirements["mpirun"],
                            "must_pass":True})

    installed.validate_stack(to_validate)

    # If we got here, reconciliation software is ready to go. Now check to
    # whether mpi can really grab the number of threads requested.
    if do_reconcile:
        check_mpi_configuration(num_threads)

    # --------------------------------------------------------------------------
    # Final sanity checks

    # If we're doing a reconciliation, make sure we can actually get placement
    # of all species on the tree.
    if do_reconcile:
        _ = topiary.df_to_species_tree(df,strict=True)

    # --------------------------------------------------------------------------
    # Deal with output directory

    # If no output directory is specified, make up a name
    if out_dir is None:
        if restart:
            err = "To use restart, you must specify an out_dir.\n"
            raise ValueError(err)
        rand = "".join([random.choice(string.ascii_letters) for _ in range(10)])
        out_dir = f"alignment_to_ancestors_{rand}"

    # Make output directory
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    else:

        # See if it is overwritable
        cannot_proceed = True
        if os.path.isdir(out_dir):
            if overwrite:
                shutil.rmtree(out_dir)
                os.mkdir(out_dir)
                cannot_proceed = False

        # If restart, we're okay having this directory around
        if restart:
            cannot_proceed = False

        # Throw error if not overwritable (overwrite = False or not dir)
        if cannot_proceed:
            err = f"\nout_dir '{out_dir}' exists. To proceed, delete or set\n"
            err += "overwrite = True.\n\n"
            raise FileExistsError(err)

    # Go into output directory
    current_dir = os.getcwd()
    os.chdir(out_dir)

    # This will count step we're on
    counter = 0

    # Find best phylogenetic model
    output = f"{counter:02d}_find-model"

    run_calc = _check_restart(output,restart)
    if run_calc:
        topiary.find_best_model(df,
                                gene_tree=starting_tree,
                                model_matrices=model_matrices,
                                model_rates=model_rates,
                                model_freqs=model_freqs,
                                model_invariant=model_invariant,
                                calc_dir=output,
                                num_threads=num_threads,
                                raxml_binary=raxml_binary)
    counter += 1

    # Generate the maximum likelihood tree without bootstraps
    prev_calculation = output
    output = f"{counter:02d}_ml-tree"

    run_calc = _check_restart(output,restart)
    if run_calc:
        topiary.generate_ml_tree(prev_calculation=prev_calculation,
                                 calc_dir=output,
                                 num_threads=num_threads,
                                 raxml_binary=raxml_binary,
                                 bootstrap=False)
    counter += 1

    # If reconciling...
    if do_reconcile:

        # Reconcile without bootstrap.
        prev_calculation = output
        output = f"{counter:02d}_reconciliation"
        run_calc = _check_restart(output,restart)
        if run_calc:
            topiary.reconcile(prev_calculation=prev_calculation,
                              calc_dir=output,
                              allow_horizontal_transfer=horizontal_transfer,
                              generax_binary=generax_binary,
                              num_threads=num_threads,
                              bootstrap=False)
        counter += 1

    # Generate ancestors
    prev_calculation = output
    output = f"{counter:02d}_ancestors"
    run_calc = _check_restart(output,restart)
    if run_calc:
        topiary.generate_ancestors(prev_calculation=prev_calculation,
                                   calc_dir=output,
                                   num_threads=num_threads,
                                   alt_cutoff=alt_cutoff)
    counter += 1

    # If we're doing bootstrap...
    if do_bootstrap:

        # Generate bootstrap replicates for the tree
        prev_calculation = output
        output = f"{counter:02d}_bootstraps"
        run_calc = _check_restart(output,restart)
        if run_calc:
            topiary.generate_bootstraps(prev_calculation=prev_calculation,
                                        calc_dir=output,
                                        num_threads=num_threads,
                                        raxml_binary=raxml_binary)
        counter += 1

    os.chdir(current_dir)

    # Create an html report for the calculation
    create_report(calculation_directory=out_dir,
                  output_directory=os.path.join(out_dir,"results"),
                  overwrite=True)
