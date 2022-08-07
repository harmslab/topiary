"""
Generate bootstrap replicates for an existing tree and then calculate boostrap
supports.
"""

import topiary

from ._raxml import run_raxml, RAXML_BINARY
from topiary._private.interface import prep_calc, write_run_information

import os, shutil, glob

def generate_bootstraps(previous_dir=None,
                        df=None,
                        model=None,
                        tree_file=None,
                        output=None,
                        overwrite=False,
                        num_threads=-1,
                        num_bootstraps=None,
                        raxml_binary=RAXML_BINARY):
    """
    Generate bootstrap replicates for an existing tree and then calculate boostrap
    supports.

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
        tree_file in newick format. Used as starting point for calculation.
        Will override tree from `previous_dir` if specified.
    output : str, optional
        output directory. If not specified, create an output directory
        with form "generate_ancestors_randomletters"
    overwrite : bool, default=False
        whether or not to overwrite existing output
    num_bootstraps : int, optional
        how many bootstrap replicates to generate. If None, use autoMRE to
        automatically infer the number of replicates given the data.
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

    # Copy files in, write out alignment, move into working directory, etc.
    result = prep_calc(previous_dir=previous_dir,
                       df=df,
                       model=model,
                       tree_file=tree_file,
                       output=output,
                       overwrite=overwrite,
                       output_base="generate_bootstraps")

    df = result["df"]
    csv_file = result["csv_file"]
    model = result["model"]
    tree_file = result["tree_file"]
    alignment_file = result["alignment_file"]
    starting_dir = result["starting_dir"]
    existing_trees = result["existing_trees"]
    start_time = result["start_time"]

    other_args = []

    # If we're doing bootstrapping
    algorithm = "--bootstrap"
    if num_bootstraps is None:
        num_bootstraps = "autoMRE"
    else:
        num_bootstraps = f"{num_bootstraps:d}"

    other_args.extend(["--bs-trees",num_bootstraps,"--bs-write-msa"])

    # Run raxml to generate bootstrap replicates
    cmd1 = run_raxml(algorithm=algorithm,
                     alignment_file=alignment_file,
                     tree_file=tree_file,
                     model=model,
                     dir_name="working",
                     seed=True,
                     num_threads=num_threads,
                     raxml_binary=raxml_binary,
                     other_args=other_args)

    # Get bs_file
    bs_file = os.path.join("working","alignment.phy.raxml.bootstraps")

    os.path.join("working/*.bootstraps")
    cmd2 = run_raxml(algorithm="--support",
                     tree_file=tree_file,
                     num_threads=1,
                     log_to_stdout=False,
                     suppress_output=True,
                     dir_name="combine-bootstraps",
                     other_files=[bs_file],
                     other_args=["--bs-trees","alignment.phy.raxml.bootstraps","--redo"])


    os.mkdir("output")

    # Copy trees from previous calculation in. This will preserve any that our
    # new calculation did not wipe out.
    for t in existing_trees:
        tree_filename = os.path.split(t)[-1]
        shutil.copy(t,os.path.join("output",tree_filename))

    # Grab the ML tree and store as tree.newick
    shutil.copy(os.path.join("input","tree.newick"),
                os.path.join("output","tree.newick"))

    # Grab tree with bootstraps and store as tree_supports.newick
    shutil.copy(os.path.join("combine-bootstraps","tree.newick.raxml.support"),
                os.path.join("output","tree_supports.newick"))

    # Copy bootstrap results to the output directory
    bs_out = os.path.join("output","bootstrap_replicates")
    os.mkdir(bs_out)
    bsmsa = glob.glob(os.path.join("working","alignment.phy.raxml.bootstrapMSA.*.phy"))
    for b in bsmsa:
        number = int(b.split(".")[-2])
        shutil.copy(b,os.path.join(bs_out,f"bsmsa_{number:04d}.phy"))
    shutil.copy(os.path.join("combine-bootstraps","alignment.phy.raxml.bootstraps"),
                os.path.join("output","bootstrap_replicates","bootstraps.newick"))

    # Write run information
    write_run_information(outdir="output",
                          df=df,
                          calc_type="ml_bootstrap",
                          model=model,
                          cmd=f"{cmd1}; {cmd2}",
                          start_time=start_time)


    print(f"\nWrote results to {os.path.abspath('output')}\n")

    # Leave working directory
    os.chdir(starting_dir)

    # Create plot holding tree
    return topiary.draw.tree(run_dir=output,
                             output_file=os.path.join(output,
                                                      "output",
                                                      "summary-tree.pdf"))