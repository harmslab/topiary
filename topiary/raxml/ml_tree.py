__description__ = \
"""
Generate maximum likelihood ancestors using raxml.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2021-07-22"

import topiary

from ._raxml import prep_calc, run_raxml, RAXML_BINARY

import os, shutil



def generate_ml_tree(previous_dir=None,
                     output=None,
                     df=None,
                     model=None,
                     tree_file=None,
                     threads=1,
                     raxml_binary=RAXML_BINARY,
                     bootstrap=True):
    """
    Generate maximum likelihood tree from an alignment given a substitution
    model.


    previous_dir: directory containing previous calculation. prep_calc will
                  grab the the csv, model, and tree from the previous run.
    output: output directory. If not specified, create an output directory with
            form "generate_ml_tree_randomletters"

    df: topiary data frame or csv written out from topiary df
    model: model (e.g. LG+G8).
    tree_file: tree_file in newick format. If not specified, a parsimony tree
               will be generated. used as starting point.

    threads: number of threads to use
    raxml_binary: what raxml binary to use
    bootstrap: whether or not to do bootstrap replicates
    """

    # Copy files in, write out alignment, move into working directory, etc.


    result = prep_calc(previous_dir=previous_dir,
                       output=output,
                       df=df,
                       model=model,
                       tree_file=tree_file,
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
        other_args.extend(["--bs-trees","autoMRE","--bs-write-msa"])

    # Run raxml to create tree
    run_raxml(algorithm="--all",
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
    shutil.copy(os.path.join("working","alignment.raxml.support"),
                os.path.join(outdir,"tree.newick"))

    # Write model to a file
    f = open(os.path.join(outdir,"model.txt"),"w")
    f.write(f"{model}\n")
    f.close()

    topiary.write_dataframe(df,os.path.join(outdir,"dataframe.csv"))

    print(f"\nWrote results to {os.path.abspath(outdir)}\n")

    # Leave working directory
    os.chdir(starting_dir)
