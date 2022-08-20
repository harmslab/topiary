"""
Reconcile gene and species trees using generax with bootstrap replicates
of the gene tree and alignments.
"""

import topiary

from topiary._private import threads
from topiary._private import Supervisor

from topiary.raxml._raxml import run_raxml
from topiary.generax._generax import setup_generax
from topiary.generax._generax import run_generax
from topiary.generax._generax import GENERAX_BINARY

import ete3

from tqdm.auto import tqdm

import os
import sys
import glob
import shutil
import copy
import tarfile
import random
import string
import subprocess
import time
import multiprocessing as mp

def _check_calc_completeness(output_directory):
    """
    Check to see how far along the calculation is.
    """

    num_directories = len(glob.glob(os.path.join(output_directory,"0*")))
    num_complete = 0

    with tqdm(total=num_directories) as pbar:
        while num_complete < num_directories:
            num_complete = len(glob.glob(os.path.join(output_directory,"0*","complete*")))
            pbar.n = num_complete
            pbar.refresh()
            time.sleep(1)


def _create_bootstrap_dirs(df,
                           model,
                           gene_tree,
                           species_tree,
                           allow_horizontal_transfer,
                           seed,
                           bootstrap_directory,
                           overwrite,
                           generax_binary):
    """
    Prepare for generax bootstrap by constructing all bootstrap directories.
    Do a dummy run of generax in each one, creating a file called run_generax.sh
    that can be invoked later.

    Parameters
    ----------
    df : pandas.DataFrame or str, optional
        topiary data frame or csv written out from topiary df. Will override
        dataframe from `prev_calculation` if specified.
    model : str, optional
        model (i.e. "LG+G8"). Will override model from `prev_calculation`
        if specified.
    gene_tree : str, optional
        gene tree in newick format.
    species_tree : str, optional
        species tree in newick format.
    allow_horizontal_transfer : bool, default=True
        whether to allow horizontal transfer during reconcilation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    seed : bool,int,str
        If true, pass a randomly generated seed to raxml. If int or str, use
        that as the seed. (passed via --seed)
    bootstrap_directory : str
        directory with gene tree bootstrap replicates
    overwrite : bool, default=False
        whether or not to overwrite existing calc_dir directory
    generax_binary : str, optional
        what generax binary to use

    Returns
    -------
    calc_dirs : list
        list of directories in which to do calculations
    """

    starting_dir = os.getcwd()

    # dataframe for boostrap replicates --> only need keep == True
    bs_df = df.loc[df.keep,:]

    # Read species tree
    if species_tree is not None:
        species_tree = species_tree
    else:
        species_tree, dropped = topiary.df_to_species_tree(bs_df,strict=True)
        species_tree.resolve_polytomy()

    # Create replicate directory
    replicate_dir = os.path.abspath(os.path.join("replicates"))
    os.mkdir(replicate_dir)

    # Find and sort bootstrap alignment files
    alignment_files = glob.glob(os.path.join(bootstrap_directory,"*.phy"))
    alignment_files = [(int(os.path.split(a)[-1].split("_")[1].split(".")[0]),a)
                       for a in alignment_files]
    alignment_files.sort()
    alignment_files = [a[1] for a in alignment_files]

    # Load bootstrap trees
    tree_files = glob.glob(os.path.join(bootstrap_directory,"*.newick"))[0]
    trees = []
    with open(tree_files) as f:
        for line in f:
            trees.append(ete3.Tree(line.strip()))

    # Sanity check on number of trees vs number of alignments
    if len(trees) != len(alignment_files):
        err = "\nNumber of bootstrap trees does not match the number of bootstrap\n"
        err += "alignments.\n\n"
        raise ValueError(err)

    template_species_tree = None
    template_control_file = None
    template_link_file = None
    template_keep_mask = None

    # For every alignment
    for i in tqdm(range(len(alignment_files))):

        rep_number = f"{(i+1):05d}"

        this_df = bs_df.copy()

        out_dir = os.path.abspath(os.path.join(replicate_dir,rep_number))

        # Write to replicate working directory
        keep_mask = setup_generax(this_df,
                                  trees[i],
                                  model,
                                  out_dir,
                                  keep_mask=template_keep_mask,
                                  species_tree=template_species_tree,
                                  mapping_link_file=template_link_file,
                                  control_file=template_control_file)

        # If this is None, we have only done first iteration. Capture those
        # files to pass into next run. (Saves a ton of time to not have to
        # recalculate species tree, traverse trees for mapping, etc.)
        if template_species_tree is None:
            template_species_tree = os.path.abspath(os.path.join(out_dir,"species_tree.newick"))
            template_control_file = os.path.abspath(os.path.join(out_dir,"control.txt"))
            template_link_file = os.path.abspath(os.path.join(out_dir,"mapping.link"))
            template_keep_mask = keep_mask.copy()

        # Copy the alignment file in, overwriting what was generated by
        # setup_generax
        align_file = os.path.join(out_dir,"alignment.phy")
        shutil.copy(alignment_files[i],align_file)

        os.chdir(replicate_dir)

        # This is a dummy run -- exactly as it would be called in a standard
        # reconciliation -- but now we're writing resulting command to
        # run_generax.sh in the run directory. Send in with 1 thread because
        # this job should run on a single thread.
        cmd = run_generax(run_directory=out_dir,
                          allow_horizontal_transfer=allow_horizontal_transfer,
                          seed=seed,
                          generax_binary=generax_binary,
                          num_threads=1,
                          log_to_stdout=False,
                          suppress_output=True,
                          write_to_script="run_generax.sh")

    os.chdir(starting_dir)

    return replicate_dir

def _run_bootstrap_calculations(replicate_dir,num_threads):
    """
    Run generax in parallel using mpirun for all directories.

    Parameters
    ----------
    replicate_dir : str
        directory with bootstrap replicates for reconcilations
    num_threads : int
        number of mpi slots (passed to mpirun -np num_threads)
    """

    print("\nGenerating reconciliation bootstraps.\n",flush=True)

    # This is a status bar that we spawn on its own thread that will spew onto
    # stderr as the calculations are completed in each directory.
    status_bar = mp.Process(target=_check_calc_completeness,
                            args=(replicate_dir,))
    status_bar.start()

    # This TMPDIR insanity is to prevent MPI from choking because temporary
    # directory names can get too long for it's buffer to handle. (Really.)
    try:
        old_tmpdir = os.environ["TMPDIR"]
    except KeyError:
        old_tmpdir = None

    os.environ["TMPDIR"] = "/tmp/"

    # Get real paths of script, mpirun, python, and replicate dir. This should
    # (hopefully!) mean we don't run into problems when we pop out onto
    # who-knows-what node...
    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_to_run = os.path.join(script_dir,"_generax_mpi_worker.py")
    script_to_run = os.path.realpath(script_to_run)

    mpirun = os.path.realpath(shutil.which("mpirun"))
    python_exec = os.path.realpath(sys.executable)
    rep_dir = os.path.realpath(replicate_dir)

    # Unique run_id so the program can recognize its own claim files (as
    # opposed to claim files left over from previous runs that may have
    # crashed)
    run_id = "".join([random.choice(string.ascii_letters) for _ in range(10)])

    # Run the worker script with mpirun if more than one thread is requested
    if num_threads > 1:
        cmd = [mpirun,"-np",f"{num_threads}"]
    else:
        cmd = []

    cmd.append(python_exec)
    cmd.append(script_to_run)
    cmd.append(rep_dir)
    cmd.append(run_id)

    # Launch command
    subprocess.run(cmd)

    # If we get here, the job is done whether the status bar is or not. Wait
    # two seconds to let the status bar finish if it's going to so the log
    # reflects the current status of the run.
    time.sleep(2)
    status_bar.kill()

    # Revert TMPDIR if necessary
    if old_tmpdir is not None:
        os.environ["TMPDIR"] = old_tmpdir


def _combine_bootstrap_calculations(replicate_dir,reconciled_tree):
    """
    Combine the results from the parallel generax bootstrap runs.

    replicate_dir : str
        directory containing the replicated generax reconcilation directories
        for the parallel bootstrap reconciliation calculation
    reconciled_tree : str
        newick file holding ML reconciled tree on which the bootstrap branch
        supports will be loaded

    Returns
    -------
    converged : bool
        whether or not the autoMRE a posterior convergence test indicates
        converged bootstrap values

    Notes
    -----
    This is not thread-safe. The generax bootstrap code uses an embarrassingly
    parallel strategy: each directory is running on its on MPI slot. This
    combining code goes through those directories and assembles the resulting
    output, even if not all directory calculations are complete. It does not
    write those directories, so it should be safe there, BUT it writes its
    results without any worry about collisions with other instances of
    _combine_results. Do not run this on more than one process.
    """

    # Paths within replicates directory we care about
    gene_tree_path = os.path.join("result",
                                  "results",
                                  "reconcile",
                                  "geneTree.newick")

    reconcile_path = os.path.join("result",
                                  "reconciliations")

    # Load strings for each bootstrap replicate into tree_stack
    tree_stack = {}
    for d in os.listdir(replicate_dir):

        rep = os.path.join(replicate_dir,d)
        if not os.path.isdir(rep):
            continue

        # See if this replicate has completed its run.
        has_run = glob.glob(os.path.join(rep,"completed_run-*"))
        if len(has_run) == 0:
            continue

        # Read file into tree_stack dict
        f = open(os.path.join(rep,gene_tree_path))
        tree_stack[d] = f.read().strip()
        f.close()

    # Sort by tree name (00001, 00002, ...)
    bs_keys = list(tree_stack.keys())
    bs_keys.sort()

    # Construct a single newick file with bootstrap replicate trees loaded in
    # line-by-line
    with open("bs-trees.newick","w") as f:
        for key in bs_keys:
            f.write(f"{tree_stack[key]}\n")

    # Combine bootstrap replicates using raxml
    if os.path.isdir("combine_with_raxml"):
        shutil.rmtree("combine_with_raxml")

    # Combine bootstrap replicates into a set of supports
    cmd = run_raxml(run_directory="combine_with_raxml",
                    algorithm="--support",
                    tree_file=reconciled_tree,
                    num_threads=1,
                    log_to_stdout=False,
                    suppress_output=True,
                    other_files=["bs-trees.newick"],
                    other_args=["--bs-trees","bs-trees.newick","--redo"])

    # Grab support tree out of raxml directory
    shutil.copy(os.path.join("combine_with_raxml",
                             "tree.newick.raxml.support"),
                "reconciled-tree_supports.newick")

    # Test for convergence
    # Combine bootstrap replicates using raxml
    if os.path.isdir("test_convergence"):
        shutil.rmtree("test_convergence")

    cmd = run_raxml(run_directory="test_convergence",
                    algorithm="--bsconverge",
                    tree_file=None,
                    num_threads=1,
                    log_to_stdout=False,
                    suppress_output=True,
                    other_files=["bs-trees.newick"],
                    other_args=["--bs-trees","bs-trees.newick"])

    converged = False
    with open(os.path.join("test_convergence","bs-trees.newick.raxml.log")) as f:
        for line in f:
            if line.startswith("Bootstrapping test converged"):
                converged = True

    return converged


def reconcile_bootstrap(df,
                        model,
                        gene_tree,
                        species_tree,
                        reconciled_tree,
                        allow_horizontal_transfer,
                        seed,
                        bootstrap_directory,
                        restart,
                        overwrite,
                        supervisor,
                        num_threads,
                        generax_binary,
                        raxml_binary):
    """
    Reconcile gene and species trees using generax with bootstrap replicates
    of the gene tree and alignments.

    Parameters
    ----------
    df : pandas.DataFrame or str, optional
        topiary data frame or csv written out from topiary df. Will override
        dataframe from `prev_calculation` if specified.
    model : str, optional
        model (i.e. "LG+G8"). Will override model from `prev_calculation`
        if specified.
    gene_tree : str, ete3.Tree, dendropy.tree, optional
        gene tree file for calculation. Will override tree in `prev_calculation`.
        If this an ete3 or dendropy tree, it will be written out with leaf
        names and branch lengths; all other data will be dropped.
    species_tree : str, ete3.Tree, dendropy.tree, optional
        species tree file for calculation. Will override tree in `prev_calculation`.
        If this an ete3 or dendropy tree, it will be written out with leaf
        names; all other data will be dropped.
    reconciled_tree : str, ete3.Tree, dendropy.tree, optional
        reconciled tree file for calculation. Will override tree in `prev_calculation`.
        If this an ete3 or dendropy tree, it will be written out with leaf
        names; all other data will be dropped.
    allow_horizontal_transfer : bool, default=True
        whether to allow horizontal transfer during reconcilation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    seed : bool,int,str
        If true, pass a randomly generated seed to raxml. If int or str, use
        that as the seed. (passed via --seed)
    bootstrap: bool, default=False
        whether or not to do bootstrap replicates. if True, prev_calculation must
        point to a raxml ml_bootstrap run
    restart : str, optional
        if specified, should point to replicate_dir for restart. restart the job
        from where it stopped. incompatible with overwrite
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

    os.chdir(supervisor.working_dir)

    if restart is None:

        # Create stack of directories
        supervisor.event("Creating bootstrap directories.",
                         model=model,
                         gene_tree=gene_tree,
                         species_tree=species_tree,
                         allow_horizontal_transfer=allow_horizontal_transfer,
                         seed=seed,
                         bootstrap_directory=bootstrap_directory,
                         overwrite=overwrite,
                         generax_binary=generax_binary)

        replicate_dir = _create_bootstrap_dirs(df=df,
                                               model=model,
                                               gene_tree=gene_tree,
                                               species_tree=species_tree,
                                               allow_horizontal_transfer=allow_horizontal_transfer,
                                               seed=seed,
                                               bootstrap_directory=bootstrap_directory,
                                               overwrite=overwrite,
                                               generax_binary=generax_binary)

    else:
        replicate_dir = restart

    supervisor.event("Running bootstrap calculations.",
                     replicate_dir=replicate_dir,
                     num_threads=num_threads)
    _run_bootstrap_calculations(replicate_dir,num_threads)

    supervisor.event("Combining bootstrap calculations.",
                     replicate_dir=replicate_dir,
                     reconciled_tree=reconciled_tree)
    converged = _combine_bootstrap_calculations(replicate_dir,reconciled_tree)

    # Grab species tree before compressing replicates
    supervisor.stash(os.path.join(replicate_dir,"00001","species_tree.newick"),
                     "species-tree.newick")

    # Compress big, complicated replicates directory and delete
    print("\nCompressing replicates.\n",flush=True)
    f = tarfile.open("replicates.tar.gz","w:gz")
    f.add("replicates")
    f.close()
    shutil.rmtree("replicates")

    # Copy in tree with supports
    supervisor.stash(os.path.join(supervisor.working_dir,"reconciled-tree_supports.newick"))

    # Record whether bootstrapping converged
    supervisor.update("bootstrap_converged",converged)

    # Write message indicating where to look for further output
    msg = "For more information on the reconcilation events (orthgroups,\n"
    msg += "event counts, full nhx files, etc.) please check the maximum\n"
    msg += "likelihood reconciliation output directory that was used as\n"
    msg += "input for this bootstrap calculation.\n"

    f = open(os.path.join(supervisor.output_dir,"reconciliations.txt"),"w")
    f.write(msg)
    f.close()

    os.chdir(supervisor.starting_dir)
    return supervisor.finalize(successful=True,plot_if_success=True)
