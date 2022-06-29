"""
Generate maximum likelihood tree from an alignment given an evolutionary
model.
"""

import topiary

from ._raxml import run_raxml, RAXML_BINARY
from topiary.external._interface import prep_calc, write_run_information

import os, shutil, glob

def generate_ml_tree(previous_dir=None,
                     df=None,
                     model=None,
                     tree_file=None,
                     output=None,
                     overwrite=False,
                     num_threads=-1,
                     raxml_binary=RAXML_BINARY,
                     bootstrap=False):
    """
    Generate maximum likelihood tree from an alignment given an evolutionary
    model.

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
    num_threads : int, default=-1
        number of threads to use. if -1, use all avaialable
    raxml_binary : str, optional
        what raxml binary to use
    bootstrap : bool, default=False
        whether or not to do bootstrap replicates

    Returns
    -------
    Python.core.display.Image or None
        if running in jupyter notebook, return Image of maximum likelihood tree;
        otherwise, return None
    """

    # Copy files in, write out alignment, move into working directory, etc.


    result = prep_calc(previous_dir=previous_dir,
                       df=df,
                       model=model,
                       tree_file=tree_file,
                       output=output,
                       overwrite=overwrite,
                       output_base="generate_ml_tree")

    df = result["df"]
    csv_file = result["csv_file"]
    model = result["model"]
    tree_file = result["tree_file"]
    alignment_file = result["alignment_file"]
    starting_dir = result["starting_dir"]

    other_args = []

    # If we're doing bootstrapping
    if bootstrap:
        algorithm = "--all"
        other_args.extend(["--bs-trees","autoMRE","--bs-write-msa"])
    else:
        algorithm = "--search"

    # Run raxml to create tree
    cmd = run_raxml(algorithm=algorithm,
                    alignment_file=alignment_file,
                    tree_file=tree_file,
                    model=model,
                    dir_name="working",
                    seed=True,
                    num_threads=num_threads,
                    raxml_binary=raxml_binary,
                    other_args=other_args)

    outdir = "output"
    os.mkdir(outdir)

    # Grab the final tree and store as tree.newick
    if bootstrap:
        shutil.copy(os.path.join("working","alignment.phy.raxml.support"),
                    os.path.join(outdir,"tree.newick"))
    else:
        shutil.copy(os.path.join("working","alignment.phy.raxml.bestTree"),
                    os.path.join(outdir,"tree.newick"))

    # Write run information
    write_run_information(outdir=outdir,
                          df=df,
                          calc_type="ml_tree",
                          model=model,
                          cmd=cmd)

    # Copy bootstrap results to the output directory
    if bootstrap:
        bs_out = os.path.join(outdir,"bootstrap_replicates")
        os.mkdir(bs_out)
        bsmsa = glob.glob(os.path.join("working","alignment.phy.raxml.bootstrapMSA.*.phy"))
        for b in bsmsa:
            number = int(b.split(".")[-2])
            shutil.copy(b,os.path.join(bs_out,f"bsmsa_{number:04d}.phy"))
        shutil.copy(os.path.join("working","alignment.phy.raxml.bootstraps"),
                    os.path.join(outdir,"bootstrap_replicates","bootstraps.newick"))

    print(f"\nWrote results to {os.path.abspath(outdir)}\n")

    # Leave working directory
    os.chdir(starting_dir)

    # Create plot holding tree
    ret = topiary.draw.ml_tree(run_dir=output,
                               output_file=os.path.join(output,
                                                        "output",
                                                        "summary-tree.pdf"))
    if topiary._in_notebook:
        return ret
