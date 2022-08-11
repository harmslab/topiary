"""
Reconcile gene and species trees using generax with bootstrap replicates
of the gene tree and alignments.
"""

import topiary

from topiary._private import threads
from topiary._private.interface import prep_calc
from topiary._private.interface import write_run_information
from topiary.raxml._raxml import run_raxml
from topiary.generax._generax import setup_generax, run_generax, GENERAX_BINARY

import ete3

from tqdm.auto import tqdm

import os
import glob
import shutil
import copy
import tarfile
import random
import string
import subprocess

def _create_bootstrap_dirs(previous_dir=None,
                           df=None,
                           model=None,
                           tree_file=None,
                           species_tree_file=None,
                           output=None,
                           overwrite=False,
                           allow_horizontal_transfer=True,
                           generax_binary=GENERAX_BINARY):
    """
    Prepare for generax bootstrap by constructing all bootstrap directories.
    Do a dummy run of generax in each one, creating a file called run_generax.sh
    that can be invoked later.

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
    species_tree_file : str
        species tree in newick format
    output: str, optional
        output directory. If not specified, create an output directory with
        form "generax_reconcilation_randomletters"
    overwrite : bool, default=False
        whether or not to overwrite existing output directory
    allow_horizontal_transfer : bool, default=True
        whether to allow horizontal transfer during reconcilation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    generax_binary : str, optional
        what generax binary to use

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
                       tree_file=tree_file,
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
    if calc_type not in ["ml_bootstrap"]:
        err = f"\nPrevious dir calc_type is '{calc_type}' but should be \n"
        err += "ml_bootstrap. Bootstrap reconciliation requires a\n"
        err +  "maximum-likelihood tree calculation directory with bootstraps\n"
        err += "as its input.\n\n"
        raise ValueError(err)

    # Make sure bootstrap directory exists
    bs_dir = os.path.join(previous_dir,"output","bootstrap_replicates")
    if not os.path.isdir(bs_dir):
        err = f"\ninput directory '{previous_dir}' does not have an\n"
        err += "output/bootstrap_replicates directory. Was this calculation run\n"
        err += "with bootstrap=True?\n\n"
        raise FileNotFoundError(err)

    # dataframe for boostrap replicates --> only need keep == True
    bs_df = df.loc[df.keep,:]

    # Read species tree
    if species_tree_file is not None:
        species_tree = species_tree_file
    else:
        species_tree, dropped = topiary.df_to_species_tree(bs_df,strict=True)
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
                                  result["model"],
                                  out_dir,
                                  keep_mask=template_keep_mask,
                                  species_tree_file=template_species_tree,
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
                          generax_binary=generax_binary,
                          num_threads=1,
                          log_to_stdout=False,
                          suppress_output=True,
                          write_to_script="run_generax.sh")

    os.chdir(starting_dir)

    return replicate_dir, result


def _combine_results(prep_output):
    """
    Combine the results from the parallel generax bootstrap runs.

    Parameters
    ----------
    prep_output : dict
        dictionary output from prep_calc (called by _prepare_for_bootstrap).
    """

    ml_tree = prep_output["tree_file"]

    # Get base output directory for the calculation
    base_dir = os.path.abspath(prep_output["output"])
    base_rep = os.path.join(base_dir,"replicates")

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

    # Combine bootstrap replicates in combine_dir using raxmls

    os.chdir(combine_dir)

    cmd = run_raxml(algorithm="--support",
                    tree_file=ml_tree,
                    num_threads=1,
                    log_to_stdout=False,
                    suppress_output=True,
                    dir_name="combine_with_raxml",
                    other_files=["bs-trees.newick"],
                    other_args=["--bs-trees","bs-trees.newick","--redo"])

    os.chdir(base_dir)

    # output directory
    os.mkdir("output")

    # Copy trees from previous calculation in. This will preserve any that our
    # new calculation did not wipe out.
    for t in prep_output["existing_trees"]:
        tree_filename = os.path.split(t)[-1]
        shutil.copy(t,os.path.join("output",tree_filename))

    # Copy in ml-tree, no supports
    shutil.copy(ml_tree,
                os.path.join("output","tree.newick"))

    # Copy ml-tree with supports into output directory
    shutil.copy(os.path.join(combine_dir,
                             "combine_with_raxml",
                             "tree.newick.raxml.support"),
                os.path.join("output","tree_supports.newick"))

    # Write message indicating where to look for further output
    msg = "For more information on the reconcilation events (orthgroups,\n"
    msg += "event counts, full nhx files, etc.) please see the maximum\n"
    msg += "likelihood reconciliation output directory that was used as\n"
    msg += "input for this bootstrap calcultion.\n"

    f = open(os.path.join("output","reconciliations.txt"),"w")
    f.write(msg)
    f.close()

    # Write run information
    write_run_information(outdir="output",
                          df=prep_output["df"],
                          calc_type="reconciliation_bootstrap",
                          model=prep_output["model"],
                          cmd=None)

    # Compress big, complicated replicates directory and delete
    print("\nCompressing replicates.\n",flush=True)
    f = tarfile.open("replicates.tar.gz","w:gz")
    f.add("replicates")
    f.close()
    shutil.rmtree("replicates")

    print(f"\nWrote results to {os.path.abspath('output')}\n",flush=True)

    os.chdir(prep_output["starting_dir"])

    # Write out a summary tree.
    return topiary.draw.tree(run_dir=base_dir,
                             output_file=os.path.join(base_dir,
                                                      "output",
                                                      "summary-tree.pdf"))

def reconcile_bootstrap(previous_dir,
                        df,
                        model,
                        tree_file,
                        species_tree_file,
                        allow_horizontal_transfer,
                        output,
                        overwrite,
                        num_threads,
                        generax_binary):
    """
    Reconcile gene and species trees using generax with bootstrap replicates
    of the gene tree and alignments. Should be invoked by
    generax.reconcile.reconcile.

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
    species_tree_file : str
        species tree in newick format
    allow_horizontal_transfer : bool, default=True
        whether to allow horizontal transfer during reconcilation. If True, use
        the "UndatedDTL" model. If False, use the "UndatedDL" model.
    output: str, optional
        output directory. If not specified, create an output directory with
        form "generax_reconcilation_randomletters"
    overwrite : bool, default=False
        whether or not to overwrite existing output directory
    num_threads : int, default=1
        number of threads to use.
    generax_binary : str, optional
        what generax binary to use

    Returns
    -------
    plot : toyplot.canvas or None
        if running in jupyter notebook, return toyplot.canvas; otherwise, return
        None.
    """

    print("Creating reconciliation bootstrap directories.\n",flush=True)

    # Create stack of directories
    replicate_dir, prep_output = _create_bootstrap_dirs(previous_dir=previous_dir,
                                                        df=df,
                                                        model=model,
                                                        tree_file=tree_file,
                                                        species_tree_file=species_tree_file,
                                                        output=output,
                                                        overwrite=overwrite)

    print("\nGenerating reconciliation bootstraps.\n",flush=True)

    run_id = "".join([random.choice(string.ascii_letters) for _ in range(10)])

    # This TMPDIR insanity is to prevent MPI from choking because temporary
    # directory names can get too long for it's buffer to handle. (Really.)
    try:
        old_tmpdir = os.environ["TMPDIR"]
    except KeyError:
        old_tmpdir = None

    os.environ["TMPDIR"] = "/tmp/"

    # Get full path of script to run
    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_to_run = os.path.join(script_dir,"_generax_mpi_worker.py")

    # Run the worker script with mpirun
    cmd = ["mpirun","-np",f"{num_threads}","python",script_to_run]
    cmd.append(replicate_dir)
    cmd.append(run_id)

    # Launch command
    subprocess.run(cmd)

    # Revert TMPDIR if necessary
    if old_tmpdir is not None:
        os.environ["TMPDIR"] = old_tmpdir

    return _combine_results(prep_output)
