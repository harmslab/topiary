"""
Pipeline that starts from an alignment, finds the best phylogenetic model,
builds a maximum likelihood tree, reconciles this tree with the species tree,
and then infers ancestral proteins.
"""

import topiary
from topiary.external.raxml import RAXML_BINARY
from topiary.external.generax import GENERAX_BINARY
from topiary._private import installed, software_requirements, check

import os, random, string, shutil

def _check_restart(output,restart):

    run_calc = True
    if restart:

        # See if json file is there. If so, assume done.
        json_file = os.path.join(output,"output","run_parameters.json")
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
                           no_reconcile=False,
                           allow_horizontal_transfer=False,
                           alt_cutoff=0.25,
                           model_matrices=["cpREV","Dayhoff","DCMut","DEN","Blosum62",
                                           "FLU","HIVb","HIVw","JTT","JTT-DCMut","LG",
                                           "mtART","mtMAM","mtREV","mtZOA","PMB",
                                           "rtREV","stmtREV","VT","WAG","LG4M","LG4X"],
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
    ancestral protein sequences.

    Parameters
    ----------
    df : pandas.DataFrame or str
        topiary data frame or csv written out from topiary df.
    out_dir : str, optional
        output directory. If not specified, create an output directory with the
        format "alignment_to_ancestors_randomletters"
    starting_tree : str, optional
        tree in newick format. This will be used for the best model
        inference and starting tree for the maximum likelihood tree estimation.
        If not specified, the maximum parsimony tree is generated and used.
    no_bootstrap : bool, default=False
        do not do bootstrap replicates
    no_reconcile : bool, default=False
        do not reconcile gene and species trees
    allow_horizontal_transfer : bool, default=False
        whether to allow horizontal transfer during reconcilation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    alt_cutoff : float, default=0.25
        cutoff to use for altAll alternate ancestral protein sequence
        generation. Should be between 0 and 1.
    model_matrices : list, default=["cpREV","Dayhoff","DCMut","DEN","Blosum62","FLU","HIVb","HIVw","JTT","JTT-DCMut","LG","mtART","mtMAM","mtREV","mtZOA","PMB","rtREV","stmtREV","VT","WAG","LG4M","LG4X"]
        list of model matrices to check. If calling from command line, these
        can be specified directly (:code:`--model_matrices LG JTT ...`) or by specifying
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
    threads : int, default=-1
        number of threads to use. if -1 use all available
    raxml_binary : str, optional
        raxml binary to use
    generax_binary : str, optional
        what generax binary to use
    """

    no_bootstrap = check.check_bool(no_bootstrap,"no_bootstrap")
    no_reconcile = check.check_bool(no_reconcile,"no_reconcile")

    # Flip logic from user interface (where flags turn off bootstrap and
    # reconcilation) to more readable flags that turn on (do_bootstrap,
    # do_reconcile)
    if no_bootstrap:
        do_bootstrap = False
    else:
        do_bootstrap = True

    if no_reconcile:
        do_reconcile = False
    else:
        do_reconcile = True
    
    overwrite = check.check_bool(overwrite,"overwrite")
    restart = check.check_bool(restart,"restart")

    if overwrite and restart:
        err = "overwrite and restart flags are incompatible.\n"
        raise ValueError(err)

    to_validate = [{"program":"raxml-ng",
                              "min_version":software_requirements["raxml-ng"],
                              "must_pass":True}]

    if do_reconcile:

        to_validate.append({"program":"generax",
                            "min_version":software_requirements["generax"],
                            "must_pass":True})
        to_validate.append({"program":"mpirun",
                            "min_version":software_requirements["mpirun"],
                            "must_pass":True})

    # Make sure the software stack is valid before doing anything
    installed.validate_stack(to_validate)

    # If no output directory is specified, make up a name
    if out_dir is None:
        if restart:
            err = "To use restart, you must specify an out_dir\n"
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

        # If restart, do some checking
        if restart:
            cannot_proceed = False

        # Throw error if not overwritable (overwrite = False or not dir)
        if cannot_proceed:
            err = f"\nout_dir '{out_dir}' exists. To proceed, delete or set\n"
            err += "overwrite = True.\n\n"
            raise FileExistsError(err)

    # If a dataframe was specified as a string, copy it in to output directory
    if issubclass(type(df),str):
        if not os.path.exists(df):
            err = f"\ndataframe '{df}' does not exist.\n\n"
            raise FileNotFoundError(err)

        df_base = os.path.split(df)[-1]
        out_df = os.path.join(out_dir,df_base)
        if not os.path.isfile(out_df):
            shutil.copy(df,out_df)
        df = df_base

    # If tree is specified as a string, copy it in to output directory
    if starting_tree is not None:
        if issubclass(type(starting_tree),str):
            if not os.path.exists(starting_tree):
                err = f"\starting_tree '{starting_tree}' does not exist.\n\n"
                raise FileNotFoundError(err)

            tree_base = os.path.split(df)[-1]
            out_tree = os.path.join(out_dir,tree_base)
            if not os.path.isfile(out_tree):
                shutil.copy(starting_tree,out_tree)
            starting_tree = tree_base

    # Go into output directory
    current_dir = os.getcwd()
    os.chdir(out_dir)


    counter = 0

    # Find best phylogenetic model
    output = f"{counter:02d}_find-model"

    run_calc = _check_restart(output,restart)
    if run_calc:
        topiary.find_best_model(df,
                                tree_file=starting_tree,
                                model_matrices=model_matrices,
                                model_rates=model_rates,
                                model_freqs=model_freqs,
                                model_invariant=model_invariant,
                                output=output,
                                num_threads=num_threads,
                                raxml_binary=raxml_binary)
    counter += 1

    # Generate the maximum likelihood tree without bootstraps
    previous_dir = output
    output = f"{counter:02d}_ml-tree"

    run_calc = _check_restart(output,restart)
    if run_calc:
        topiary.generate_ml_tree(previous_dir=previous_dir,
                                 output=output,
                                 num_threads=num_threads,
                                 raxml_binary=raxml_binary,
                                 bootstrap=False)
    counter += 1

    # If reconciling...
    if do_reconcile:

        # Reconcile without bootstrap.
        previous_dir = output
        output = f"{counter:02d}_reconciliation"
        run_calc = _check_restart(output,restart)
        if run_calc:
            topiary.reconcile(previous_dir=previous_dir,
                              output=output,
                              allow_horizontal_transfer=allow_horizontal_transfer,
                              generax_binary=generax_binary,
                              num_threads=num_threads,
                              use_mpi=False, #### HACK HACK HACK
                              bootstrap=False)
        counter += 1

    # Generate ancestors
    previous_dir = output
    output = f"{counter:02d}_ancestors"
    run_calc = _check_restart(output,restart)
    if run_calc:
        topiary.generate_ancestors(previous_dir=previous_dir,
                                   output=output,
                                   num_threads=num_threads,
                                   alt_cutoff=alt_cutoff)
    counter += 1

    # If we're doing bootstrap, reconcile all bootstrap replicates and
    # generate ancestor with branch supports
    if do_bootstrap:

        # Generate bootstrap replicates for the tree
        previous_dir = "01_ml-tree"
        output = f"{counter:02d}_bootstraps"
        run_calc = _check_restart(output,restart)
        if run_calc:
            topiary.generate_bootstraps(previous_dir=previous_dir,
                                        output=output,
                                        num_threads=num_threads,
                                        raxml_binary=raxml_binary)
        counter += 1

        # Generate bootstraps for reconciliation
        if do_reconcile:
            previous_dir = output
            output = f"{counter:02d}_reconciliation-bootstraps"
            run_calc = _check_restart(output,restart)
            if run_calc:
                topiary.reconcile(previous_dir=previous_dir,
                                  output=output,
                                  allow_horizontal_transfer=allow_horizontal_transfer,
                                  generax_binary=generax_binary,
                                  num_threads=num_threads,
                                  use_mpi=False,
                                  bootstrap=do_bootstrap)
            counter += 1

        # Generate final ancestors on tree with branch supports
        previous_dir = output
        output = f"{counter:02d}_ancestors_with_branch_supports"
        run_calc = _check_restart(output,restart)
        if run_calc:
            topiary.generate_ancestors(previous_dir=previous_dir,
                                       output=output,
                                       num_threads=num_threads,
                                       alt_cutoff=alt_cutoff)


    os.chdir(current_dir)
