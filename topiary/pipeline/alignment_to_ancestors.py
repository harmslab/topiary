"""
Pipeline that starts from an alignment, finds the best phylogenetic model,
builds a maximum likelihood tree, reconciles this tree with the species tree,
and then infers ancestral proteins.
"""

import topiary
from topiary.external.raxml import RAXML_BINARY
from topiary.external.generax import GENERAX_BINARY

import os, random, string, shutil

def alignment_to_ancestors(df,
                           out_dir=None,
                           starting_tree=None,
                           do_bootstrap=True,
                           allow_horizontal_transfer=False,
                           alt_cutoff=0.25,
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
    do_bootstrap : bool, default=True
        whether or not to do bootstrap replicates
    allow_horizontal_transfer : bool, default=False
        whether to allow horizontal transfer during reconcilation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    alt_cutoff : float, default=0.25
        cutoff to use for altAll alternate ancestral protein sequence
        generation. Should be between 0 and 1.
    overwrite : bool, default=False
        whether or not to overwrite existing output
    threads : int, default=-1
        number of threads to use. if -1 use all available
    raxml_binary : str, optional
        raxml binary to use
    generax_binary : str, optional
        what generax binary to use


    ## NOT SURE WHETHER/HOW TO INCLUDE...
    model_matrices : list, optional
        list of model matrices to check
    model_rates : list, optional
        ways to treat model rates
    model_freqs : list, optional
        ways to treat model freqs.
    model_invariant : list, optional
        ways to treat invariant alignment columns

    Returns
    -------
    None
    """

    # If no output directory is specified, make up a name
    if out_dir is None:
        rand = "".join([random.choice(string.ascii_letters) for _ in range(10)])
        out_dir = f"alignment_to_ancestors_{rand}"

    # If output directory already exists
    if os.path.exists(out_dir):

        # See if it is overwritable
        cannot_overwrite = True
        if os.path.isdir(out_dir):
            if overwrite:
                shutil.rmtree(out_dir)
                cannot_overwrite = False

        # Throw error if not overwritable (overwrite = False or not dir)
        if cannot_overwrite:
            err = f"\nout_dir '{out_dir}' exists. To proceed, delete or set\n"
            err += "overwrite = True.\n\n"
            raise FileExistsError(err)

    # Make output directory
    os.mkdir(out_dir)

    # If a dataframe was specified as a string, copy it in to output directory
    if issubclass(type(df),str):
        if not os.path.exists(df):
            err = f"\ndataframe '{df}' does not exist.\n\n"
            raise FileNotFoundError(err)

        df_base = os.path.split(df)[-1]
        out_df = os.path.join(out_dir,df_base)
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
            shutil.copy(starting_tree,out_tree)
            starting_tree = tree_base


    # Go into output directory
    current_dir = os.getcwd()
    os.chdir(out_dir)

    topiary.find_best_model(df,
                            tree_file=starting_tree,
                            # model_matrices
                            # model_rates
                            # model_freqs
                            # model_invariant
                            output="00_find-model",
                            threads=num_threads,
                            raxml_binary=raxml_binary)

    topiary.generate_ml_tree(previous_dir="00_find-model",
                             output="01_ml-tree",
                             threads=num_threads,
                             raxml_binary=raxml_binary,
                             bootstrap=do_bootstrap)

    topiary.reconcile(previous_dir="01_ml-tree",
                      output="02_reconciliation",
                      allow_horizontal_transfer=allow_horizontal_transfer,
                      generax_binary=generax_binary)

    topiary.generate_ancestors(previous_dir="02_reconciliation",
                               output="03_ancestors",
                               alt_cutoff=alt_cutoff)

    os.chdir(current_dir)

    ## model_matrices=["cpREV","Dayhoff","DCMut","DEN","Blosum62",
    ##                 "FLU","HIVb","HIVw","JTT","JTT-DCMut","LG",
    ##                 "mtART","mtMAM","mtREV","mtZOA","PMB",
    ##                 "rtREV","stmtREV","VT","WAG","LG4M","LG4X",
    ##                 "PROTGTR"],
    ## model_rates=["","G8"],
    ## model_freqs=["","FC","FO"],
    ## model_invariant=["","IO","IC"],
