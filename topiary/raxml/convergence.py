"""
Check a newick file containing bootstrap replicates for convergence using
a specified convergence cutoff.
"""
from topiary._private import check
from topiary._private.interface import rmtree
from ._raxml import run_raxml, RAXML_BINARY

import pandas as pd

import shutil
import os

def _parse_convergence_file(converge_file):
    """
    Parse the output of a --bsconverge RAxML-NG call.

    Parameters
    ----------
    converge_file : str
        convergence log file

    Returns
    -------
    converged : bool
        whether or not the bootstraps are converged
    df : pandas.DataFrame
        dataframe holding convergence results
    """

    out = {"trees":[],
           "avg_wrf":[],
           "avg_wrf_pct":[],
           "perms_below_cutoff":[],
           "converged":[]}

    collecting = False
    with open(converge_file) as f:
        for line in f:
            if not collecting and line.startswith(" # trees "):
                collecting = True
                continue

            if collecting and line.startswith("Bootstopping"):
                break

            if collecting:
                col = line.split()
                out["trees"].append(int(col[0].strip()))
                out["avg_wrf"].append(float(col[1].strip()))
                out["avg_wrf_pct"].append(float(col[2].strip()))
                out["perms_below_cutoff"].append(int(col[3].strip()))
                if col[4].strip() == "NO":
                    out["converged"] = False
                else:
                    out["converged"] = True

    df = pd.DataFrame(out)
    converged = df.loc[:,"converged"].iloc[-1]

    return converged, df


def check_convergence(bs_newick,
                      converge_cutoff=0.03,
                      seed=True,
                      calc_dir=None,
                      num_threads=1,
                      raxml_binary=RAXML_BINARY):
    """
    Check a newick file containing bootstrap replicates for convergence using
    a specified convergence cutoff.

    bs_newick : str
        newick file containing bootstrap replicate trees
    converge_cutoff : float, default=0.03
        convergence cutoff. permutations with WRF below cutoff that stop
        bootstrapping
    seed : bool,int,str
        If true, pass a randomly generated seed to raxml. If int or str, use
        that as the seed. (passed via --seed)
    calc_dir : str, optional
        write output to this directory. if None, create temporary directory,
        parse output, and then delete temporary directory
    num_threads : int, default=1
        run calculation on num_threads threads
    raxml_binary : str, optional
        what raxml binary to use

    Returns
    -------
    converged : bool
        whether or not the bootstraps are converged
    df : pandas.DataFrame
        dataframe holding convergence results
    """

    # Process/check input arguments

    bs_newick = f"{bs_newick}"
    if not os.path.isfile(bs_newick):
        err = f"\nbs_newick '{bs_newick}' does not exist.\n\n"
        raise FileNotFoundError(err)
    bs_newick_basename = os.path.basename(bs_newick)

    converge_cutoff = check.check_float(converge_cutoff,
                                        "converge_cutoff",
                                        minimum_allowed=0,
                                        maximum_allowed=1)

    num_threads = check.check_int(num_threads,
                                  "num_threads",
                                  minimum_allowed=1)

    delete_output = False
    if calc_dir is None:
        calc_dir = "tmp_test_conv"
        delete_output = True
    calc_dir = f"{calc_dir}"


    # Run convergence test
    cmd = run_raxml(run_directory=calc_dir,
                    algorithm="--bsconverge",
                    tree_file=None,
                    seed=seed,
                    num_threads=num_threads,
                    log_to_stdout=False,
                    suppress_output=True,
                    other_files=[bs_newick],
                    other_args=["--bs-trees",bs_newick_basename,
                                "--bs-cutoff",f"{converge_cutoff}"])

    # Parse results
    conv_file = os.path.join(calc_dir,f"{bs_newick_basename}.raxml.log")
    converged, df = _parse_convergence_file(conv_file)

    # Clean up
    if delete_output:
        rmtree(calc_dir)

    return converged, df
