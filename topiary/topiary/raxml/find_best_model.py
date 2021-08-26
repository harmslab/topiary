__description__ = \
"""
Find ML model using parsimony tree.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2021-07-22"

from ._raxml import prep_calc, run_raxml, RAXML_BINARY, gen_seed

import pandas as pd
import numpy as np

import os, re, shutil

def _parse_raxml_info_for_aic(info_file):
    """
    Open a log file from a tree evaluation run and get log likelihood,
    number of fit parameters, and various AIC scores.
    """

    # Open the file and read lines
    out = {}
    with open(info_file,'r') as f:
        for line in f:

            # Look for likelihood
            if re.search("Final LogLikelihood:",line):
                out["L"] = float(line.strip().split(":")[-1])
                continue

            # Look for number of parameters
            if re.search("Free parameters",line):
                out["N"] = int(line.strip().split(":")[-1])
                continue

            if re.search("AIC score",line):
                cols = line.split("/")
                for c in cols:

                    value = float(c.strip().split(":")[-1])
                    key = c.strip().split(":")[0].split()[0].strip()
                    out[key] = value

    # Return L, N, and AIC
    return out

def _generate_parsimony_tree(alignment_file,
                             dir_name="parsimony-tree",
                             threads=1,
                             raxml_binary=RAXML_BINARY):
    """
    Generate a parsimony tree from an alignment.

    alignment_file: alignment file in .phy format
    dir_name: name to give directory
    threads: number of threads to use
    raxml_binary: raxml binary to use
    """

    run_raxml(algorithm="--start",
              alignment_file=alignment_file,
              dir_name=dir_name,
              seed=True,
              model="LG",
              threads=threads,
              raxml_binary=raxml_binary,
              log_to_stdout=False,
              other_args=["--tree","pars{1}"])

def find_best_model(df,
                    tree_file=None,
                    model_matrices=["Blosum62","cpREV"], #,"Dayhoff","DCMut","DEN",
                                    #"FLU","HIVb","HIVw","JTT","JTT-DCMut","LG",
                                    #"mtART","mtMAM","mtREV","mtZOA","PMB",
                                    #"rtREV","stmtREV","VT","WAG","LG4M","LG4X",
                                    #"PROTGTR"],
                    model_rates=["","G8"],
                    model_freqs=["","FC","FO"],
                    model_invariant=["","IO","IC"],
                    output=None,
                    threads=1,
                    raxml_binary=RAXML_BINARY):
    """
    Find the best phylogentic model to use for tree and ancestor reconstruction
    given an alignment and (possibly) a tree.

    df: topiary data frame or csv written out from topiary df
    tree_file: tree file in newick format. If not specified, parsimony tree
               is generated and used
    model_matrices: list of model matrices to check
    model_rates: ways to treat model rates
    model_freqs: ways to treat model freqs.
    output: directory for output. it none, will generate random name
    threads: number of threads to use
    raxml_binary: raxml binary to use
    """

    # Copy files in, write out alignment, move into working directory, etc.
    result = prep_calc(df=df,
                       output=output,
                       other_files=[tree_file],
                       output_base="find_best_model")

    df = result["df"]
    csv_file = result["csv_file"]
    alignment_file = result["alignment_file"]
    tree_file = result["other_files"][0]
    starting_dir = result["starting_dir"]

    # Generate a parsimony tree if not was specified
    if tree_file is None:
        _generate_parsimony_tree(alignment_file,
                                 dir_name="01_make-parsimony-tree",
                                 threads=threads,
                                 raxml_binary=raxml_binary)
        tree_file = "02_parsimony-tree.newick"
        shutil.copy(os.path.join("01_make-parsimony-tree",
                                 "alignment.raxml.startTree"),
                    tree_file)

    # Dictionary to hold stats for each model
    out = {"model":[]}

    seed = gen_seed()

    # All possible models, dropping rate, freq, invariant for LG4M and LG4X.
    num_mat = len([m for m in model_matrices if m not in ["LG4M","LG4X"]])
    num_models = num_mat*len(model_rates)*len(model_freqs)*len(model_invariant) + 2

    # Go over all combos of the requested matrices, rates, and freqs.
    model_counter = 1
    for matrix in model_matrices:
        for rate in model_rates:
            for freq in model_freqs:
                for invariant in model_invariant:

                    # Check for incompatible matrix/freq/rate combos
                    if matrix in ["LG4M","LG4X"]:
                        if rate != "" or freq != "" or invariant != "":
                            continue

                    # Build model string (for example: LG+G8+FC+IO)
                    model = [matrix,rate,freq,invariant]
                    model = [m for m in model if m != ""]
                    model = "+".join(model)

                    # Print model number we're trying
                    print(f"{model} ({model_counter}/{num_models})")
                    model_counter += 1

                    # Optimize branch lengths etc. on the existing tree
                    run_raxml(algorithm="--evaluate",
                              alignment_file=alignment_file,
                              tree_file=tree_file,
                              model=model,
                              seed=seed,
                              dir_name="tmp",
                              threads=threads,
                              raxml_binary=raxml_binary,
                              log_to_stdout=False)

                    # Grab the info file from this run
                    os.chdir("tmp")

                    # Get results from the info file
                    result = _parse_raxml_info_for_aic("alignment.raxml.log")
                    out["model"].append(model)
                    for r in result:
                        try:
                            out[r].append(result[r])
                        except KeyError:
                            out[r] = [result[r]]

                    # Get out of temporary directory and nuke
                    os.chdir("..")
                    shutil.rmtree("tmp")

    # Create a csv file sorted best to worst aicc
    final_df = pd.DataFrame(out)

    min_aic = np.min(final_df.AICc)
    final_df["p"] = np.exp((min_aic - final_df.AICc)/2)
    indexer = np.argsort(final_df.p)[::-1]
    final_df = final_df.iloc[indexer,:]
    final_df.to_csv("model-comparison.csv")

    # Get best model
    best_model = final_df.model.iloc[0]

    # Write model to a file
    f = open("best-model.txt","w")
    f.write(best_model)
    f.close()

    # Print best model to stdout
    print(f"\n\nBest model: {best_model}\nAICc Prob:{final_df.p.iloc[0]}\n\n")

    # Leave the output directory
    os.chdir(starting_dir)
