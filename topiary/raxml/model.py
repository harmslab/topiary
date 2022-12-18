"""
Find the best phylogentic model to use for tree and ancestor reconstruction
given an alignment and (possibly) a tree.
"""

import topiary

from ._raxml import RAXML_BINARY, run_raxml

from topiary._private import check
from topiary._private import threads
from topiary._private import Supervisor
from topiary._private import run_cleanly

import pandas as pd
import numpy as np

import os, re, shutil, random, string
from tqdm.auto import tqdm


def _generate_parsimony_tree(alignment_file,
                             run_directory="parsimony-tree",
                             seed=True,
                             supervisor=None,
                             num_threads=1,
                             raxml_binary=RAXML_BINARY):
    """
    Generate a parsimony tree from an alignment.

    Parameters
    ----------
        alignment_file: alignment file in .phy format
        run_directory: name to give directory
        num_threads: number of threads to use
        raxml_binary: raxml binary to use

    Return
    ------
        None
    """

    run_raxml(run_directory=run_directory,
              algorithm="--start",
              alignment_file=alignment_file,
              seed=seed,
              model="LG",
              log_to_stdout=False,
              suppress_output=True,
              other_args=["--tree","pars{1}"],
              supervisor=supervisor,
              num_threads=num_threads,
              raxml_binary=raxml_binary)

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

    # Run raxml. We catch with RuntimeError because sometimes a model/dataset
    # combo is so bad it can't optimize. Drop those from further analysis.
    current_dir = os.getcwd()
    try:
        run_raxml(**kwargs)
        success = True
    except RuntimeError:
        success = False

    # This makes sure we came back to starting directory, specifically in
    # event of a raxml crash
    os.chdir(current_dir)

    # Get results from the info file
    tmp_dir = kwargs["run_directory"]
    info_file = os.path.join(tmp_dir,"alignment.phy.raxml.log")

    if success:
        result = _parse_raxml_info_for_aic(info_file)
    else:
        result = None

    # Nuke temporary directory
    shutil.rmtree(tmp_dir)

    return result

@run_cleanly
def find_best_model(df,
                    gene_tree=None,
                    model_matrices=["cpREV","Dayhoff","DCMut","DEN","Blosum62",
                                    "FLU","HIVb","HIVw","JTT","JTT-DCMut","LG",
                                    "mtART","mtMAM","mtREV","mtZOA","PMB",
                                    "rtREV","stmtREV","VT","WAG"],
                    model_rates=["","G8"],
                    model_freqs=["","FC","FO"],
                    model_invariant=["","IC","IO"],
                    seed=None,
                    calc_dir="find_best_model",
                    overwrite=False,
                    supervisor=None,
                    num_threads=-1,
                    raxml_binary=RAXML_BINARY):
    """
    Find the best phylogentic model to use for tree and ancestor reconstruction
    given an alignment and (possibly) a tree.

    Parameters
    ----------
    df : pandas.DataFrame or str
        topiary data frame or csv written out from topiary df.
    gene_tree : str
        gene_tree in newick format. If not specified, parsimony tree is
        generated and used
    model_matrices : list, default=["cpREV","Dayhoff","DCMut","DEN","Blosum62","FLU","HIVb","HIVw","JTT","JTT-DCMut","LG","mtART","mtMAM","mtREV","mtZOA","PMB","rtREV","stmtREV","VT","WAG"]
        list of model matrices to check
    model_rates : list, default=["","G8"]
        ways to treat model rates. If None, do not include a model rate param.
    model_freqs : list, default=["","FC","FO"]
        ways to treat model freqs. If None, use the matrix frequencies.
    model_invariant : list, default=["","IC","IO"]
        ways to treat invariant alignment columns. If None, do not have an
        invariant class.
    seed : bool,int,str
        If true, pass a randomly generated seed to raxml. If int or str, use
        that as the seed. (passed via --seed)
    calc_dir : str, default="find_best_model"
        calculation directory.
    overwrite : bool, default=False
        whether or not to overwrite existing output
    supervisor : Supervisor, optional
        instance of Supervisor for managing calculation inputs and outputs
    num_threads : int, default=-1
        number of threads to use. if -1 use all available
    raxml_binary : str, optional
        raxml binary to use

    Returns
    -------
    None
    """

    # Create supervisor if needed
    if supervisor is None:
        supervisor = Supervisor(seed=seed)

    # Create directory with the supervisor
    supervisor.create_calc_dir(calc_dir=calc_dir,
                               calc_type="find_best_model",
                               overwrite=overwrite,
                               df=df,
                               gene_tree=gene_tree)

    supervisor.check_required(required_files=["alignment.phy","dataframe.csv"])

    os.chdir(supervisor.working_dir)

    # Check sanity of model matrices, rates, freqs, invariant, num_threads
    model_matrices = check.check_iter(model_matrices,
                                      "model_matrices",
                                      required_value_type=str)

    # Deal with model rates
    if model_rates is None:
        model_rates = [""]

    model_rates = check.check_iter(model_rates,
                                   "model_rates",
                                   required_value_type=str)
    model_rates = list(model_rates)
    if "" not in model_rates:
        model_rates.insert(0,"")

    # Deal with model freqs
    if model_freqs is None:
        model_freqs = [""]

    model_freqs = check.check_iter(model_freqs,
                                   "model_freqs",
                                   required_value_type=str)
    model_freqs = list(model_freqs)
    if "" not in model_freqs:
        model_freqs.insert(0,"")

    # Deal with model invariant
    if model_invariant is None:
        model_invariant = [""]

    model_invariant = check.check_iter(model_invariant,
                                       "model_invariant",
                                       required_value_type=str)
    model_invariant = list(model_invariant)
    if "" not in model_invariant:
        model_invariant.insert(0,"")

    num_threads = threads.get_num_threads(num_threads)

    # If we already have a tree, use it
    if supervisor.gene_tree is not None:
        gene_tree = supervisor.gene_tree

    # Generate a parsimony tree if not was specified
    else:

        print("\nGenerating maximum parsimony tree.\n",flush=True)

        _generate_parsimony_tree(supervisor.alignment,
                                 run_directory="max-parsimony",
                                 seed=supervisor.seed,
                                 supervisor=supervisor,
                                 num_threads=num_threads,
                                 raxml_binary=raxml_binary)
        gene_tree = os.path.join("max-parsimony",
                                 "alignment.phy.raxml.startTree")

    print("Constructing set of possible models.")

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

                    # Make temporary directory
                    rand = "".join([random.choice(string.ascii_letters)
                                    for _ in range(10)])
                    run_directory = f"tmp_{rand}"

                    # set up kwargs for raxml calcs. single thread per calc
                    kwargs = {"run_directory":run_directory,
                              "algorithm":"--evaluate",
                              "alignment_file":supervisor.alignment,
                              "tree_file":gene_tree,
                              "model":model,
                              "seed":supervisor.seed,
                              "num_threads":1,
                              "raxml_binary":raxml_binary,
                              "log_to_stdout":False,
                              "suppress_output":True}
                    kwargs_list.append({"kwargs":kwargs})
                    models.append(model)


    supervisor.event(f"Testing {len(kwargs_list)} models.\n")

    out_list = threads.thread_manager(kwargs_list,
                                      _model_thread_function,
                                      num_threads)

    # Go through output list and store results in out
    out = {"model":[]}
    for i, result in enumerate(out_list):

        # Failed calculation
        if result is None:
            continue

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
    final_df = final_df.loc[final_df.index[indexer],:]

    # Write dataframe comparing models
    final_df.to_csv(os.path.join(supervisor.output_dir,
                                 "model-comparison.csv"))

    # Print best model to stdout
    print("\nTop 10 models:\n")
    print(f"{'model':>20s}{'AICc prob':>20s}",flush=True)
    for i in range(10):
        if i >= len(final_df):
            break
        print(f"{final_df.loc[final_df.index[i],'model']:>20s}{final_df.loc[final_df.index[i],'p']:>20.3f}",flush=True)

    best_model = final_df.loc[final_df.index[0],"model"]

    # Record model in supervisor so it will appear in run_parameters.json
    supervisor.update("model",best_model)

    # Close out
    os.chdir(supervisor.starting_dir)
    return supervisor.finalize(successful=True,plot_if_success=False)
