"""
Reconcile gene and species trees using generax with bootstrap replicates
of the gene tree and alignments.
"""

import topiary

from topiary._private import threads
from topiary.external._interface import prep_calc, write_run_information
from topiary.external.raxml._raxml import run_raxml
from topiary.external.generax._generax import setup_generax, run_generax, GENERAX_BINARY

import ete3

from tqdm.auto import tqdm

import os
import glob
import shutil
import copy

def _prepare_for_bootstrap(previous_dir=None,
                           df=None,
                           model=None,
                           tree_file=None,
                           output=None,
                           overwrite=False):
    """
    Prepare for generax bootstrap by constructing all bootstrap directories and
    then returning those directories as a list.

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
    output: str, optional
        output directory. If not specified, create an output directory with
        form "generax_reconcilation_randomletters"
    overwrite : bool, default=False
        whether or not to overwrite existing output directory

    Returns
    -------
    calc_dirs : list
        list of directories in which to do calculations
    """

    previous_dir = os.path.abspath(previous_dir)

    # Prepare for the calculation, loading in previous calculation and
    # combining with arguments as passed in.
    result = prep_calc(previous_dir=previous_dir,
                       df=df,
                       model=model,
                       tree_file=model,
                       output=output,
                       overwrite=overwrite,
                       output_base="generax_bootstrap_reconcilation")

    df = result["df"]
    csv_file = result["csv_file"]
    model = result["model"]
    tree_file = result["tree_file"]
    calc_type = result["calc_type"]
    alignment_file = result["alignment_file"]
    starting_dir = result["starting_dir"]

    # Make sure previous run has required files
    required = [df,model,tree_file]
    for r in required:
        if r is None:
            err = "\nA dataframe, model, and tree are required for this "
            err += "calculation.\n\n"
            raise ValueError(err)

    # Read previous run directory
    if calc_type != "ml_tree":
        err = f"\nPrevious dir calc_type is '{calc_type}' but should be 'ml_tree'.\n"
        err += "Bootstrap reconciliation must have a maximum-likelihood tree\n"
        err += "calculation directory as its input.\n\n"
        raise ValueError(err)

    # Make sure bootstrap directory exists
    bs_dir = os.path.join(previous_dir,"output","bootstrap_replicates")
    if not os.path.isdir(bs_dir):
        err = f"\ninput directory '{previous_dir}' does not have an\n"
        err += "output/bootstrap_replicates directory. Was this calculation run\n"
        err += "with bootstrap=True?\n\n"
        raise FileNotFoundError(err)

    bs_df = df.loc[df.keep,:]

    # Load species tree from opentree
    species_tree = topiary.get_species_tree(bs_df)
    species_tree.resolve_polytomy()

    # Create replicate directory
    replicate_dir = os.path.abspath(os.path.join("replicates"))
    os.mkdir(replicate_dir)

    # Find and sort bootstrap alignment files
    alignment_files = glob.glob(os.path.join(bs_dir,"*.phy"))
    alignment_files = [(int(os.path.split(a)[-1].split("_")[1].split(".")[0]),a)
                       for a in alignment_files]
    alignment_files.sort()
    alignment_files = [a[1] for a in alignment_files]

    # Load bootstrap trees
    tree_file = glob.glob(os.path.join(bs_dir,"*.newick"))[0]
    trees = []
    with open(tree_file) as f:
        for line in f:
            trees.append(ete3.Tree(line.strip()))

    # Sanity check on number of trees vs number of alignments
    if len(trees) != len(alignment_files):
        err = "\nNumber of bootstrap trees does not match the number of bootstrap\n"
        err += "alignments.\n\n"
        raise ValueError(err)

    # For every alignment
    calc_dirs = []
    for i in tqdm(range(len(alignment_files))):

        rep_number = f"{(i+1):05d}"

        this_df = bs_df.copy()

        out_dir = os.path.abspath(os.path.join(replicate_dir,rep_number))
        calc_dirs.append(os.path.join(out_dir))

        # Write to replicate working directory
        setup_generax(this_df,
                      trees[i],
                      result["model"],
                      out_dir,
                      species_tree=species_tree)
        os.chdir(out_dir)

        tree_file = os.path.abspath(os.path.join("tree.newick"))
        trees[i].write(format=6,outfile=tree_file)

        align_file = os.path.join("alignment.phy")
        shutil.copy(alignment_files[i],align_file)

        # Write a temporary fasta file to read back into dataframe alignment
        # column.
        first_line = True
        fasta_lines = []
        with open(alignment_files[i]) as f:

            for line in f:

                # Skip first line
                if first_line:
                    first_line = False
                    continue

                # Skip blank lines
                if line.strip() == "":
                    continue

                # Read uid and column (assumes xxxxxxxxxx seq ... format)
                col = line.strip().split()
                name = col[0]
                seq = col[1]

                fasta_lines.append(f">{name}")
                fasta_lines.append(seq)

        # Load fasta alignment into the dataframe
        this_df = topiary.io.read_fasta_into(this_df,
                                             fasta_lines,
                                             load_into_column="alignment")

        # Write dataframe
        df_file = "dataframe.csv"
        topiary.write_dataframe(this_df,df_file)

        os.chdir(replicate_dir)

    # Create directory for maximum likelihood reconcilation
    os.chdir(starting_dir)
    ml_dir = os.path.abspath(os.path.join(replicate_dir,"ml"))
    calc_dirs.append(ml_dir)

    setup_generax(df,tree_file,model,ml_dir,species_tree=species_tree)

    return result, calc_dirs

def _construct_args(calc_dirs,
                    allow_horizontal_transfer=True,
                    seed=None,
                    generax_binary=GENERAX_BINARY,
                    other_args=[],
                    num_threads=-1,
                    test_num_cores=None):
    """
    Construct arguments to pass to each thread in the pool.

    Parameters
    ----------
    calc_dirs : list
        list of directories in which to do calculation
    allow_horizontal_transfer : bool, default=False
        whether or not to allow horizontal gene transfer. This corresponds to
        the UndatedDTL (horizontal) vs UndatedDL (no horizontal) models
    seed : bool or int or str, optional
        If true, pass a randomly generated seed to generax. If int or str, use
        that as the seed (passed via --seed).
    generax_binary : str, optional
        generax binary to use
    other_args : list, optional
        other arguments to pass to generax
    num_threads : int, default=-1
        number of threads. if -1, use all available
    test_num_cores : int, optional
        for the number of cores to be test_num_cores (for testing)

    Returns
    ------
        list of args to pass for each calculation, number of threads
    """

    num_threads = threads.get_num_threads(num_threads,test_num_cores)

    kwargs_list = []
    for c in calc_dirs:
        kwargs_list.append({"run_directory":c,
                            "allow_horizontal_transfer":allow_horizontal_transfer,
                            "seed":seed,
                            "generax_binary":generax_binary,
                            "other_args":copy.deepcopy(other_args)})

    return kwargs_list, num_threads


def _generax_thread(run_directory,
                    allow_horizontal_transfer,
                    seed,
                    generax_binary,
                    other_args):
    """
    Run a generax calculation on a single thread.

    Parameters
    ----------
    run_directory : str
        directory in which to do calculation
    allow_horizontal_transfer : bool, default=False
        whether or not to allow horizontal gene transfer. This corresponds to
        the UndatedDTL (horizontal) vs UndatedDL (no horizontal) models
    seed : bool or int or str, optional
        If true, pass a randomly generated seed to generax. If int or str, use
        that as the seed (passed via --seed).
    generax_binary : str, optional
        generax binary to use
    log_to_stdout : bool, default=True
        capture log and write to std out.
    other_args : list, optional
        other arguments to pass to generax
    """

    cmd = run_generax(run_directory=run_directory,
                      allow_horizontal_transfer=allow_horizontal_transfer,
                      seed=seed,
                      generax_binary=generax_binary,
                      log_to_stdout=False,
                      other_args=other_args)

def _combine_results(prep_output):
    """
    Combine the results from the parallel generax bootstrap runs.

    Parameters
    ----------
    prep_output : dict
        dictionary output from prep_calc (called by _prepare_for_bootstrap).
    """

    # Get base output directory for the calculation
    base_dir = prep_output["output"]
    base_rep = os.path.join(base_dir,"replicates")

    # output directory
    outdir = os.path.join(base_dir,"output")
    os.mkdir(outdir)

    # combine-bootstraps directory
    combine_dir = os.path.join(base_dir,"combine-bootstraps")
    os.mkdir(combine_dir)

    # Paths within replicates directory we care about
    gene_tree_path = os.path.join("result",
                                  "results",
                                  "reconcile",
                                  "geneTree.newick")

    reconcile_path = os.path.join("result",
                                  "reconciliations")

    # Load strings for each bootstrap replicate into tree_stack
    tree_stack = {}
    for d in os.listdir(base_rep):

        # skip ML replicate
        if d == "ml":
            continue

        # Read file into tree_stack dict
        f = open(os.path.join(base_rep,d,gene_tree_path))
        tree_stack[d] = f.read().strip()
        f.close()

    # Sort by tree name (00001, 00002, ...)
    bs_keys = list(tree_stack.keys())
    bs_keys.sort()

    # Construct a single newick file with bootstrap replicate trees loaded in
    # line-by-line
    with open(os.path.join(combine_dir,"bs-trees.newick"),"w") as f:
        for key in bs_keys:
            f.write(f"{tree_stack[key]}\n")

    # Copy ml tree into combine dir
    shutil.copy(os.path.join(base_rep,"ml",gene_tree_path),
                os.path.join(combine_dir,"ml-tree.newick"))

    # Combine bootstrap replicates in combine_dir using raxmls
    starting_dir = os.getcwd()
    os.chdir(combine_dir)

    cmd = run_raxml(algorithm="--support",
                    tree_file="ml-tree.newick",
                    num_threads=1,
                    log_to_stdout=False,
                    dir_name="combine_with_raxml",
                    other_files=["bs-trees.newick"],
                    other_args=["--bs-trees","bs-trees.newick","--redo"])
    os.chdir(starting_dir)

    # Copy ml-tree with supports into output directory
    shutil.copy(os.path.join(combine_dir,
                             "combine_with_raxml",
                             "tree.newick.raxml.support"),
                os.path.join(outdir,"tree.newick"))

    # Copy reconcilation directory into output directory
    shutil.copytree(os.path.join(base_rep,"ml",reconcile_path),
                    os.path.join(outdir,"reconcilations"))

    # Get outgroups (e.g. leaves descending from each half after root)
    reconcile_file = os.path.join(base_rep,
                                  "ml",reconcile_path,
                                  "reconcile_events.newick")
    reconcile_tree = ete3.Tree(reconcile_file,format=1)
    root = reconcile_tree.get_tree_root()
    root_children = root.get_children()
    outgroup = [[n.name for n in r.get_leaves()] for r in root_children]

    # Write run information
    write_run_information(outdir=outdir,
                          df=prep_output["df"],
                          calc_type="reconciliation_bootstrap",
                          model=prep_output["model"],
                          cmd=None,
                          outgroup=outgroup)

    print(f"\nWrote results to {os.path.abspath(outdir)}\n")

    # Write out a summary tree.
    return topiary.draw.reconciliation_tree(run_dir=base_dir,
                                            output_file=os.path.join(base_dir,
                                                                     "output",
                                                                     "summary-tree.pdf"))

def reconcile_bootstrap(previous_dir=None,
                        df=None,
                        model=None,
                        tree_file=None,
                        allow_horizontal_transfer=True,
                        output=None,
                        overwrite=False,
                        num_threads=-1,
                        generax_binary=GENERAX_BINARY):
    """
    Reconcile gene and species trees using generax with bootstrap replicates
    of the gene tree and alignments.

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

    print("Creating bootstrap directories.",flush=True)

    # Create stack of directories
    prep_output, calc_dirs = _prepare_for_bootstrap(previous_dir=previous_dir,
                                                    df=df,
                                                    model=model,
                                                    tree_file=tree_file,
                                                    output=output,
                                                    overwrite=overwrite)

    # Construct keyword arguments to pass to thread
    kwargs_list, num_threads = _construct_args(calc_dirs,num_threads)

    # Run multi-threaded generax--one reconciliation per thread
    threads.thread_manager(kwargs_list,
                           _generax_thread,
                           num_threads,
                           progress_bar=True,
                           pass_lock=False)

    return _combine_results(prep_output)