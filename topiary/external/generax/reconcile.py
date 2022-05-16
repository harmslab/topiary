
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

    # Load gene tree
    gene_tree = topiary.util.load_tree(tree_file)

    # Create generax data structures
    gene_tree, species_tree, link_dict = _generax._create_generax_input(df,gene_tree)

    # Write out generax input
    _generax._write_generax_input(df,gene_tree,species_tree,link_dict,model)

    _generax.run_generax(".",
                         allow_horizontal_transfer=allow_horizontal_transfer,
                         generax_binary=generax_binary)


    # outdir = "output"
    # os.mkdir(outdir)
    #
    # # Write out a pretty version of the tree
    # shutil.copy(os.path.join("working","alignment.raxml.support"),
    #             os.path.join(outdir,"tree.newick"))
    #
    # # Write model to a file
    # f = open(os.path.join(outdir,"model.txt"),"w")
    # f.write(f"{model}\n")
    # f.close()
    #
    # topiary.write_dataframe(df,os.path.join(outdir,"dataframe.csv"))
    #
    # # Copy bootstrap results to the output directory
    # if bootstrap:
    #     bs_out = os.path.join(outdir,"bootstrap_replicates")
    #     os.mkdir(bs_out)
    #     bsmsa = glob.glob(os.path.join("working","alignment.raxml.bootstrapMSA.*.phy"))
    #     for b in bsmsa:
    #         number = int(b.split(".")[-2])
    #         shutil.copy(b,os.path.join(bs_out,f"bsmsa_{number:04d}.phy"))
    #     shutil.copy(os.path.join("working","alignment.raxml.bootstraps"),
    #                 os.path.join(outdir,"bootstrap_replicates","bootstraps.newick"))
    #
    # print(f"\nWrote results to {os.path.abspath(outdir)}\n")

    # Leave working directory
    os.chdir(starting_dir)
