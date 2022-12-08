"""
Reconcile gene and species trees using generax with bootstrap replicates
of the gene tree and alignments.
"""

import topiary

from topiary._private import Supervisor
from topiary._private import threads
from topiary._private import mpi

from topiary.raxml._raxml import run_raxml
from topiary.raxml import check_convergence
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
import pathlib
import multiprocessing as mp


def _progress_bar(replicate_dir):
    """
    Check to see how far along the calculation is.

    Parameters
    ----------
    replicate_dir : str
        directory containing replicates
    """

    total_calcs = len(glob.glob(os.path.join(replicate_dir,"0*")))

    num_complete = 0
    skipped_added = False
    with tqdm(total=total_calcs) as pbar:
        while num_complete < total_calcs:

            num_complete = len(glob.glob(os.path.join(replicate_dir,"0*","completed")))
            num_running = len(glob.glob(os.path.join(replicate_dir,"0*","running")))
            num_skipped = len(glob.glob(os.path.join(replicate_dir,"0*","skipped")))

            # Num skipped just added in
            if num_skipped > 0 and not skipped_added:

                skipped_added = True
                total_calcs = num_complete + num_running

                pbar.reset(total_calcs)
                pbar.update(num_complete)
                pbar.refresh()

            pbar.n = num_complete
            pbar.refresh()
            time.sleep(0.5)

def _check_convergence(replicate_dir,
                       converge_cutoff,
                       lock=None):
    """
    Check for convergence and indicate to all threads to terminate if true.

    Parameters
    ----------
    replicate_dir : str
        directory containing replicates
    converge_cutoff : float
        bootstrap convergence criterion. passed to --bs-cutoff
    lock : multiprocessing.Manager().Lock(), optional
        lock to control file

    Returns
    -------
    converged : bool
        whether or not the bootstraps are converged
    df : pandas.DataFrame
        dataframe holding convergence results
    """

    if lock is None:
        lock = threads.MockLock()

    # Grab copy of the newick file with lock to avoid collisons with threads
    # that might be writing to it.
    lock.acquire()
    try:
        shutil.copy(os.path.join(replicate_dir,"bs-trees.newick"),"tmp.newick")
    finally:
        lock.release()

    # Don't do calculation of the number of trees is only one
    f = open('tmp.newick')
    lines = [line for line in f.readlines() if line.strip() != ""]
    f.close()

    if len(lines) < 2:
        os.remove("tmp.newick")
        return False, None

    # Check for convergence
    converged, df = check_convergence("tmp.newick",
                                      converge_cutoff=converge_cutoff)

    # Delete temporary newick file
    os.remove("tmp.newick")

    # If calculation has converged...
    if converged:

        # Get list of directories
        dirs = [os.path.join(replicate_dir,d) for d in os.listdir(replicate_dir)]
        dirs = [d for d in dirs if os.path.isdir(d)]
        dirs.sort()

        num_completed = 0
        num_running = 0

        # Lock to prevent other threads from starting a calculation or updating
        # the progress bar while this is running.
        lock.acquire()
        try:

            # Go through all directories
            for d in dirs:

                # Is calculation done?
                if os.path.isfile(os.path.join(d,"completed")):
                    num_completed += 1

                # if not...
                else:

                    # If calculation running?
                    if os.path.isfile(os.path.join(d,"running")):
                        num_running += 1

                    # if not, write file telling other threads to skip
                    else:
                        pathlib.Path(os.path.join(d,"skipped")).touch()

        # Release lock
        finally:
            lock.release()

    return converged, df


def _generax_thread_function(replicate_dir,
                             converge_cutoff,
                             is_manager,
                             hosts,
                             lock=None):
    """
    Run a generax calculation in parallel, checking for and avoiding collisions
    with other workers.

    Parameters
    ----------

    replicate_dir : str
        directory containing replicates
    converge_cutoff : float
        bootstrap convergence criterion. passed to --bs-cutoff
    is_manager : bool
        whether or not this is the manager thread that will check for
        convergence.
    hosts : list
        list of hosts on which to run the calculation. passed to mpirun via
        --hosts ",".join(hosts)
    lock : multiprocessing.Manager().Lock()
        lock to allow multiple threads to access files
    """

    if lock is None:
        lock = threads.MockLock()

    # Change into replicate_directory and get sorted list of bootstrap
    # replicates.
    os.chdir(replicate_dir)
    dirs = [d for d in os.listdir(".") if os.path.isdir(d)]
    dirs.sort()

    # Construct a base mpirun command that the generax commands will be
    # appended to
    base_cmd = ["mpirun","--host",",".join(hosts)]

    # Path to result tree within each directory
    result_tree = os.path.join("result","results","reconcile","geneTree.newick")

    # Go through each directory
    converged = None
    df = None
    for d in dirs:

        # figure out what directory to go into. Use the lock to make sure only
        # one process is staking a claim at a time.
        lock.acquire()
        try:

            # If this directory has already been run
            if os.path.isfile(os.path.join(d,"completed")):
                continue

            # If this directory is already running
            if os.path.isfile(os.path.join(d,"running")):
                continue

            # If this directory has already been set to skip
            if os.path.isfile(os.path.join(d,"skipped")):
                continue

            # Stake a claim
            os.chdir(d)
            pathlib.Path("running").touch()

        finally:
            lock.release()

        # If we got all the way here, this directory is ours to run in.

        # Read run_generax.sh script and construct an mpirun command with the
        # contents of the scripts
        f = open("run_generax.sh")
        contents = f.read()
        f.close()

        # Split on "\n" and strip redirect if presents
        bash_cmd = contents.split("\n")[0].strip()
        bash_cmd = bash_cmd.split("&>")[0]
        bash_cmd = bash_cmd.split()
        bash_cmd = [c.strip() for c in bash_cmd if c.strip() != ""]
        cmd = base_cmd[:]
        cmd.extend(bash_cmd)

        # Launch as a subprocess
        ret = subprocess.run(cmd,capture_output=True)

        # Write stdout and stderr
        f = open("stdout.log","w")
        f.write(ret.stdout.decode())
        f.close()

        f = open("stderr.log","w")
        f.write(ret.stderr.decode())
        f.close()

        # If failure, raise a RuntimeError
        if ret.returncode != 0:
            err = f"\ngenerax crashed in directory {d}. Writing stderr and\n"
            err += "stout there.\n\n"
            raise RuntimeError(err)

        # Grab result tree
        f = open(result_tree,'r')
        tree = f.read().strip()
        f.close()

        # Update status bar and final tree file
        lock.acquire()
        try:

            # Update replicates tree file
            f = open(os.path.join("..","bs-trees.newick"),"a")
            f.write(f"{tree}\n")
            f.close()

        finally:
            lock.release()

        # Create a completed file, nuke running file, move to next directory.
        pathlib.Path("completed").touch()
        os.remove("running")
        os.chdir("..")

        # For the manager thread, check for convergence
        if is_manager:
            converged, df = _check_convergence(replicate_dir=".",
                                               converge_cutoff=converge_cutoff,
                                               lock=lock)
            if converged:
                break

    os.chdir("..")

    if is_manager:
        return (converged, df)

    return None



def _build_replicate_dirs(df,
                          model,
                          gene_tree,
                          species_tree,
                          allow_horizontal_transfer,
                          seed,
                          bootstrap_directory,
                          overwrite,
                          generax_binary):
    """
    Prepare for generax bootstrap by constructing individual directories that
    hold each bootstrap replicate. Do a dummy run of generax in each one,
    creating a file called run_generax.sh that can be invoked later.

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
        whether to allow horizontal transfer during reconciliation. If True, use
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
    replicate_dir : str
        directory containing replicate directories
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
        # run_generax.sh in the run directory.
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

def _clean_replicate_dir(replicate_dir):
    """
    Remove previous run information.

    Parameters
    ----------
    replicate_dir : str
        directory containing replicates
    """

    running = glob.glob(os.path.join(replicate_dir,"0*","running"))
    for f in running:
        os.remove(f)

    skipped = glob.glob(os.path.join(replicate_dir,"0*","skipped"))
    for f in skipped:
        os.remove(f)


def _construct_args(replicate_dir,
                    converge_cutoff,
                    num_threads,
                    threads_per_rep):
    """
    Construct a list of arguments to pass to each thread in the pool.

    Parameters
    ----------
    replicate_dir: str
        directory containing replicates
    converge_cutoff : float
        bootstrap convergence criterion. passed to --bs-cutoff
    num_threads : int
        total number of mpi slots to use for the calculation
    threads_per_rep : int
        number of slots to use per replicate.

    Returns
    -------
    kwargs_list : list
        list of dictionaries with kwargs to pass for each calculation
    num_threads : int
        number of total calculations to start via thread_manager
    """

    hosts = mpi.get_hosts(num_threads)

    kwargs_list = []
    for i in range(0,len(hosts),threads_per_rep):

        if i == 0:
            is_manager = True
        else:
            is_manager = False

        kwargs_list.append({"replicate_dir":replicate_dir,
                            "is_manager":is_manager,
                            "hosts":hosts[i:i+threads_per_rep],
                            "converge_cutoff":converge_cutoff})

    return kwargs_list, len(kwargs_list)


def _run_bootstrap_calculations(replicate_dir,
                                converge_cutoff,
                                num_threads,
                                threads_per_rep):
    """
    Run generax in parallel using mpirun for all directories.

    Parameters
    ----------
    replicate_dir : str
        directory with bootstrap replicates for reconciliations
    converge_cutoff : float
        bootstrap convergence criterion. passed to --bs-cutoff
    num_threads : int
        number of parallel jobs to start
    threads_per_rep : int
        number of threads to use per replicate. only used if bootstrap = True
    """

    print("\nGenerating reconciliation bootstraps.\n",flush=True)

    # This is a status bar that we spawn on its own thread that will spew onto
    # stderr as the calculations are completed in each directory.
    status_bar = mp.Process(target=_progress_bar,args=(replicate_dir,))
    status_bar.start()

    # This TMPDIR insanity is to prevent MPI from choking because temporary
    # directory names can get too long for it's buffer to handle. (Really.)
    try:
        old_tmpdir = os.environ["TMPDIR"]
    except KeyError:
        old_tmpdir = None
    os.environ["TMPDIR"] = "/tmp/"


    kwargs_list, num_threads = _construct_args(replicate_dir,
                                               converge_cutoff,
                                               num_threads,
                                               threads_per_rep)

    # Launch calculation.
    try:
        results = threads.thread_manager(kwargs_list,
                                         _generax_thread_function,
                                         num_threads=num_threads,
                                         progress_bar=False,
                                         pass_lock=True)

    # This clean up step makes sure we kill the status bar thread if the
    # calculation crashes.
    except Exception as e:
        status_bar.kill()
        raise e

    # If we get here, the job is done whether the status bar is or not. Wait
    # two seconds to let the status bar finish if it's going to so the log
    # reflects the current status of the run.
    time.sleep(2)
    status_bar.kill()

    # Revert TMPDIR if necessary
    if old_tmpdir is not None:
        os.environ["TMPDIR"] = old_tmpdir

    # Only the first thread will spit out results: a tuple with convergence
    # status (True or False) and the dataframe corresponding to last convergence
    # test. Get those results. If df is None, convergence test not run yet;
    # run it.
    converged, df = results[0]
    if df is None:
        converged, df = _check_convergence(replicate_dir=replicate_dir,
                                           converge_cutoff=converge_cutoff)

    return converged, df


def reconcile_bootstrap(df,
                        model,
                        gene_tree,
                        species_tree,
                        reconciled_tree,
                        allow_horizontal_transfer,
                        seed,
                        bootstrap_directory,
                        converge_cutoff,
                        restart,
                        overwrite,
                        supervisor,
                        num_threads,
                        threads_per_rep,
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
        whether to allow horizontal transfer during reconciliation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    seed : bool,int,str
        If true, pass a randomly generated seed to raxml. If int or str, use
        that as the seed. (passed via --seed)
    bootstrap: bool, default=False
        whether or not to do bootstrap replicates. if True, prev_calculation must
        point to a raxml ml_bootstrap run
    converge_cutoff : float, default=0.03
        bootstrap convergence criterion. This is RAxML-NG default, passed to
        --bs-cutoff.
    restart : str, optional
        if specified, should point to replicate_dir for restart. restart the job
        from where it stopped. incompatible with overwrite
    overwrite : bool, default=False
        whether or not to overwrite existing calc_dir directory
    supervisor : Supervisor, optional
        supervisor instance to keep track of inputs and outputs
    num_threads : int, default=-1
        number of threads to use. if -1 use all available.
    threads_per_rep : int, default=1
        number of threads to use per replicate.
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

        replicate_dir = _build_replicate_dirs(df=df,
                                              model=model,
                                              gene_tree=gene_tree,
                                              species_tree=species_tree,
                                              allow_horizontal_transfer=allow_horizontal_transfer,
                                              seed=seed,
                                              bootstrap_directory=bootstrap_directory,
                                              overwrite=overwrite,
                                              generax_binary=generax_binary)

    else:

        # Restart, making sure any leftover claim/skipped files are removed
        replicate_dir = restart
        _clean_replicate_dir(replicate_dir)



    supervisor.event("Running bootstrap calculations.",
                     replicate_dir=replicate_dir,
                     converge_cutoff=converge_cutoff,
                     num_threads=num_threads,
                     threads_per_rep=threads_per_rep)
    converged, df = _run_bootstrap_calculations(replicate_dir,
                                                converge_cutoff,
                                                num_threads,
                                                threads_per_rep)

    # Write convergence report and whether this converged or not
    df.to_csv("bootstrap-convergence-report.csv")
    supervisor.stash("bootstrap-convergence-report.csv")
    supervisor.update("bootstrap_converged",bool(converged))

    # Combine bootstrap replicates into a set of supports
    supervisor.event("Combining bootstrap calculations.",
                     replicate_dir=replicate_dir,
                     reconciled_tree=reconciled_tree)

    bs_trees = os.path.join(replicate_dir,"bs-trees.newick")
    cmd = run_raxml(run_directory="combine_with_raxml",
                    algorithm="--support",
                    tree_file=reconciled_tree,
                    num_threads=1,
                    log_to_stdout=False,
                    suppress_output=True,
                    other_files=[bs_trees],
                    other_args=["--bs-trees","bs-trees.newick","--redo"])

    # Copy reconciled tree with suppports into output
    supervisor.stash(os.path.join("combine_with_raxml","tree.newick.raxml.support"),
                     target_name="reconciled-tree_supports.newick")

    # Grab species tree before compressing replicates
    supervisor.stash(os.path.join(replicate_dir,"00001","species_tree.newick"),
                     "species-tree.newick")

    # Compress big, complicated replicates directory and delete
    print("\nCompressing replicates.\n",flush=True)
    f = tarfile.open("replicates.tar.gz","w:gz")
    f.add("replicates")
    f.close()
    shutil.rmtree("replicates")

    # Write message indicating where to look for further output
    msg = "For more information on the reconciliation events (orthgroups,\n"
    msg += "event counts, full nhx files, etc.) please check the maximum\n"
    msg += "likelihood reconciliation output directory that was used as\n"
    msg += "input for this bootstrap calculation.\n"

    f = open(os.path.join(supervisor.output_dir,"reconciliations.txt"),"w")
    f.write(msg)
    f.close()

    os.chdir(supervisor.starting_dir)
    return supervisor.finalize(successful=True,plot_if_success=True)
