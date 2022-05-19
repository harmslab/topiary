__description__ = \
"""
Generate maximum likelihood ancestors using raxml.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2021-07-22"

import topiary

from ._raxml import run_raxml, RAXML_BINARY
from topiary.external.interface import prep_calc

import os, shutil, glob

def generate_ml_tree(previous_dir=None,
                     df=None,
                     model=None,
                     tree_file=None,
                     output=None,
                     overwrite=False,
                     threads=-1,
                     raxml_binary=RAXML_BINARY,
                     bootstrap=False):
    """
    Generate maximum likelihood tree from an alignment given a substitution
    model.

    previous_dir: directory containing previous calculation. prep_calc will
                  grab the the csv, model, and tree from the previous run.
    df: topiary data frame or csv written out from topiary df
    model: model (e.g. LG+G8).
    tree_file: tree_file in newick format. If not specified, a parsimony tree
               will be generated. used as starting point.
    output: output directory. If not specified, create an output directory with
            form "generate_ml_tree_randomletters"
    overwrite: whether or not to overwrite existing output (default False)
    threads: number of threads to use. if -1 use all available
    raxml_binary: what raxml binary to use
    bootstrap: whether or not to do bootstrap replicates
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
    run_raxml(algorithm=algorithm,
              alignment_file=alignment_file,
              tree_file=tree_file,
              model=model,
              dir_name="working",
              seed=True,
              threads=threads,
              raxml_binary=raxml_binary,
              other_args=other_args)

    outdir = "output"
    os.mkdir(outdir)

    # Write out a pretty version of the tree
    if bootstrap:
        shutil.copy(os.path.join("working","alignment.phy.raxml.support"),
                    os.path.join(outdir,"tree.newick"))
    else:
        shutil.copy(os.path.join("working","alignment.phy.raxml.bestTree"),
                    os.path.join(outdir,"tree.newick"))

    # Write model to a file
    f = open(os.path.join(outdir,"model.txt"),"w")
    f.write(f"{model}\n")
    f.close()

    topiary.write_dataframe(df,os.path.join(outdir,"dataframe.csv"))

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

    try:

        df.nickname
        tip_columns = ["species","nickname"]
        local_df = None

    except AttributeError:

        local_df = df.copy()
        trunc_name = []
        for i in range(len(local_df)):
            if len(local_df["name"].iloc[i]) > 10:
                trunc_name.append(f"{local_df['name'].iloc[i][:10]}...")
            else:
                trunc_name.append(f"{local_df['name'].iloc[i]}")

        local_df["trunc_name"] = trunc_name
        tip_columns = ["species","trunc_name"]

    ret = topiary.draw.ml_tree(ml_dir=output,
                               df=local_df,
                               output=os.path.join(output,
                                                   "output",
                                                   "summary-tree.pdf"),
                               tip_columns=tip_columns)

    if topiary._in_notebook:
        return ret
