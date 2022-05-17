
import topiary
from topiary.external.interface import prep_calc

from . import _generax

import numpy as np

import os, glob, shutil

def reconcile(previous_dir=None,
              output=None,
              df=None,
              model=None,
              tree_file=None,
              allow_horizontal_transfer=False,
              generax_binary="generax"):


    result = prep_calc(previous_dir=previous_dir,
                       output=output,
                       df=df,
                       model=model,
                       tree_file=model,
                       output_base="generax_reconcilation")

    df = result["df"]
    csv_file = result["csv_file"]
    model = result["model"]
    tree_file = result["tree_file"]
    alignment_file = result["alignment_file"]
    starting_dir = result["starting_dir"]

    _generax.setup_generax(df,tree_file,model,"working")

    _generax.run_generax(run_directory="working",
                         allow_horizontal_transfer=allow_horizontal_transfer,
                         generax_binary=generax_binary)

    outdir = "output"
    os.mkdir(outdir)

    shutil.copy(os.path.join("working","result","results","reconcile","geneTree.newick"),
                os.path.join("output","tree.newick"))

    # Write model to a file
    f = open(os.path.join(outdir,"model.txt"),"w")
    f.write(f"{model}\n")
    f.close()

    # Write dataframe
    topiary.write_dataframe(df,os.path.join(outdir,"dataframe.csv"))

    # Copy reconcilation information
    shutil.copytree(os.path.join("working","result","reconciliations"),
                    os.path.join("output","reconcilations"))

    # Leave working directory
    os.chdir(starting_dir)
