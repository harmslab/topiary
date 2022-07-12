"""
Find the best phylogentic model to use for tree and ancestor reconstruction
given an alignment and (possibly) a tree.
"""

import topiary

from ._raxml import RAXML_BINARY, run_raxml
from topiary._private.interface import prep_calc, gen_seed, write_run_information
from topiary._private import check
from topiary._private import threads

import pandas as pd
import numpy as np

import os, re, shutil, sys, random, string
import multiprocessing as mp
from tqdm.auto import tqdm


def _generate_parsimony_tree(alignment_file,
                             dir_name="parsimony-tree",
                             num_threads=1,
                             raxml_binary=RAXML_BINARY):
    """
    Generate a parsimony tree from an alignment.

    Parameters
    ----------
        alignment_file: alignment file in .phy format
        dir_name: name to give directory
        num_threads: number of threads to use
        raxml_binary: raxml binary to use

    Return
    ------
        None
    """

    run_raxml(algorithm="--start",
              alignment_file=alignment_file,
              dir_name=dir_name,
              seed=True,
              model="LG",
              num_threads=num_threads,
              raxml_binary=raxml_binary,
              log_to_stdout=False,
              other_args=["--tree","pars{1}"])


def _parse_raxml_info_for_aic(info_file):
    """
    Open a log file from a tree evaluation run and get log likelihood,
    number of fit parameters, and various AIC scores.

    Parameters
    ----------
        info_file: raxml info file

    Return
    ------
        dictionary holding likelihood, number of free parameters, and aic score(s)
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


def _model_thread_function(kwargs):
    """
    Run raxml on a single thread and extract likelihood, number of parameters,
    and AIC scores.

    Parameters
    ----------
    kwargs : dict
        keyword arguments to pass to run_raxml

    Returns
    -------
    result : dict
        dictionary with "L", "N", and various AIC scores.
    """

    # Run raxml
    run_raxml(**kwargs)

    # Get results from the info file
    tmp_dir = kwargs["dir_name"]
    info_file = os.path.join(tmp_dir,"alignment.phy.raxml.log")
    result = _parse_raxml_info_for_aic(info_file)

    # Nuke temporary directory
    shutil.rmtree(tmp_dir)

    return result


def find_best_model(df,
                    tree_file=None,
                    model_matrices=["cpREV","Dayhoff","DCMut","DEN","Blosum62",
                                    "FLU","HIVb","HIVw","JTT","JTT-DCMut","LG",
                                    "mtART","mtMAM","mtREV","mtZOA","PMB",
                                    "rtREV","stmtREV","VT","WAG","LG4M","LG4X"],
                    model_rates=["","G8"],
                    model_freqs=["","FC","FO"],
                    model_invariant=["","IC","IO"],
                    output=None,
                    overwrite=False,
                    num_threads=-1,
                    raxml_binary=RAXML_BINARY):
    """
    Find the best phylogentic model to use for tree and ancestor reconstruction
    given an alignment and (possibly) a tree.

    Parameters
    ----------
    df : pandas.DataFrame or str
        topiary data frame or csv written out from topiary df.
    tree_file : str
        tree_file in newick format. If not specified, parsimony tree is
        generated and used
    model_matrices : list, default=["cpREV","Dayhoff","DCMut","DEN","Blosum62","FLU","HIVb","HIVw","JTT","JTT-DCMut","LG","mtART","mtMAM","mtREV","mtZOA","PMB","rtREV","stmtREV","VT","WAG","LG4M","LG4X"]
        list of model matrices to check
    model_rates : list, default=["","G8"]
        ways to treat model rates
    model_freqs : list, default=["","FC","FO"]
        ways to treat model freqs.
    model_invariant : list, default=["","IC","IO"]
        ways to treat invariant alignment columns
    output : str, optional
        output directory. If not specified, create an output directory with
        form "find_best_model_randomletters"
    overwrite : bool, default=False
        whether or not to overwrite existing output
    num_threads : int, default=-1
        number of threads to use. if -1 use all available
    raxml_binary : str, optional
        raxml binary to use

    Returns
    -------
    None
    """

    # Copy files in, write out alignment, move into working directory, etc.
    result = prep_calc(df=df,
                       tree_file=tree_file,
                       output=output,
                       overwrite=overwrite,
                       output_base="find_best_model")

    df = result["df"]
    csv_file = result["csv_file"]
    alignment_file = result["alignment_file"]
    tree_file = result["tree_file"]
    starting_dir = result["starting_dir"]

    # Generate a parsimony tree if not was specified
    if tree_file is None:

        print("\nGenerating maximum parsimony tree.\n",flush=True)

        _generate_parsimony_tree(alignment_file,
                                 dir_name="working",
                                 num_threads=num_threads,
                                 raxml_binary=raxml_binary)
        tree_file = os.path.join("working",
                                 "alignment.phy.raxml.startTree")

    # Deal with model matrices
    model_matrices = check.check_iter(model_matrices,
                                      "model_matrices",
                                      required_value_type=str)

    # Deal with model rates
    model_rates = check.check_iter(model_rates,
                                   "model_rates",
                                   required_value_type=str)
    model_rates = list(model_rates)
    if "" not in model_rates:
        model_rates.insert(0,"")

    # Deal with model freqs
    model_freqs = check.check_iter(model_freqs,
                                   "model_freqs",
                                   required_value_type=str)
    model_freqs = list(model_freqs)
    if "" not in model_freqs:
        model_freqs.insert(0,"")

    # Deal with model invariant
    model_invariant = check.check_iter(model_invariant,
                                       "model_invariant",
                                       required_value_type=str)
    model_invariant = list(model_invariant)
    if "" not in model_invariant:
        model_invariant.insert(0,"")

    num_threads = threads.get_num_threads(num_threads)

    seed = gen_seed()

    # Go over all combos of the requested matrices, rates, freqs, invariant and
    # create kwargs for run_raxml calls
    kwargs_list = []
    models = []
    for matrix in model_matrices:
        for rate in model_rates:
            for freq in model_freqs:
                for invariant in model_invariant:

                    # Build model string (for example: LG+G8+FC+IO)
                    model = [matrix,rate,freq,invariant]
                    model = [m for m in model if m != ""]
                    model = "+".join(model)

                    # Check for incompatible matrix/freq/rate combos
                    if matrix in ["LG4M","LG4X"]:
                        if rate != "" or freq != "" or invariant != "":
                            print(f"skipping incompatible model combination {model}",flush=True)
                            continue

                    # Make temporary directory
                    rand = "".join([random.choice(string.ascii_letters)
                                    for _ in range(10)])
                    dir_name = f"tmp_{rand}"

                    kwargs = {"algorithm":"--evaluate",
                              "alignment_file":alignment_file,
                              "tree_file":tree_file,
                              "model":model,
                              "seed":seed,
                              "dir_name":dir_name,
                              "num_threads":1, # single thread per calc
                              "raxml_binary":raxml_binary,
                              "log_to_stdout":False}
                    kwargs_list.append({"kwargs":kwargs})
                    models.append(model)


    out_list = threads.thread_manager(kwargs_list,
                                      _model_thread_function,
                                      num_threads)


    # Go through output list and store results in out
    out = {"model":[]}
    for i, result in enumerate(out_list):

        out["model"].append(models[i])
        for r in result:
            try:
                out[r].append(result[r])
            except KeyError:
                out[r] = [result[r]]

    # Create a dataframe sorted best to worst aicc
    final_df = pd.DataFrame(out)

    min_aic = np.min(final_df.AICc)
    final_df["p"] = np.exp((min_aic - final_df.AICc)/2)
    final_df["p"] = final_df["p"]/np.sum(final_df["p"])
    indexer = np.argsort(final_df.p)[::-1]
    final_df = final_df.iloc[indexer,:]

    # All output goes into the output directory
    outdir = "output"
    os.mkdir(outdir)

    # Write run information
    write_run_information(outdir=outdir,
                          df=df,
                          calc_type="find_best_model",
                          model=final_df.model.iloc[0],
                          cmd=None)

    # Write dataframe comparing models
    final_df.to_csv(os.path.join(outdir,"model-comparison.csv"))

    # Print best model to stdout
    print("\nTop 10 models:\n")
    print(f"{'model':>20s}{'AICc prob':>20s}",flush=True)
    for i in range(10):
        if i >= len(final_df):
            break
        print(f"{final_df.model.iloc[i]:>20s}{final_df.p.iloc[i]:>20.3f}",flush=True)

    print(f"\nWrote results to {os.path.abspath(outdir)}\n")

    # Leave the output directory
    os.chdir(starting_dir)
