"""
Run raxml. Creates a working directory, copies in the relevant files, runs
there, and then returns to the previous directory.
"""

# raxml binary to use it not specified by user
RAXML_BINARY = "raxml-ng"

import topiary
from topiary._private import interface
from topiary._private import threads

import pandas as pd
import multiprocessing as mp
import os

def run_raxml(algorithm=None,
              alignment_file=None,
              tree_file=None,
              model=None,
              dir_name=None,
              seed=None,
              num_threads=-1,
              raxml_binary=RAXML_BINARY,
              log_to_stdout=True,
              other_args=None,
              other_files=None):
    """
    Run raxml. Creates a working directory, copies in the relevant files, runs
    there, and then returns to the previous directory.

    Parameters
    ----------
    algorithm : str
        algorithm to run (--all, --ancestral, etc.)
    alignment_file : str
        alignment file in .phy format (passed via --msa)
    tree_file : str
        tree file in .newick format (passed via --tree)
    model : str
        model in format recognized by --model
    dir_name : str,optional
        If specified, this will be the name of the working directory.
    seed : bool,int,str
        If true, pass a randomly generated seed to raxml. If int or str, use
        that as the seed. (passed via --seed)
    num_threads : int, default=-1
        number of threads (passed via --threads). if -1, use all available.
    raxml_binary : str, default=RAXML_BINARY
        raxml binary to use
    log_to_stdout : book, default=True
        capture log and write to std out.
    other_args : list-like, optional
        list of arguments to pass to raxml
    other_files : list-like, optional
        list of files to copy into working directory (besides tree_file and
        alignment_file)

    Return
    ------
    raxml_command : string
        command passed to raxml-ng as a string
    """

    # Create directory in which to do calculation
    dir_name = interface.create_new_dir(dir_name=dir_name)

    # Copy alignment and tree files into the directory (if specified)
    if alignment_file is not None:
        alignment_file = interface.copy_input_file(alignment_file,
                                                   dir_name,
                                                   file_name="alignment.phy",
                                                   put_in_input_dir=False)
    if tree_file is not None:
        tree_file = interface.copy_input_file(tree_file,
                                              dir_name,
                                              file_name="tree.newick",
                                              put_in_input_dir=False)

    # Copy in any other required files, if requested
    if other_files is not None:
        for i in range(len(other_files)):
            tail = os.path.split(other_files[i])[-1]
            other_files[i] = interface.copy_input_file(other_files[i],
                                                       dir_name,
                                                       file_name=tail,
                                                       put_in_input_dir=False)

    # Build a command list
    cmd = [raxml_binary]

    if algorithm is not None:
        cmd.append(algorithm)

    if alignment_file is not None:
        cmd.extend(["--msa",alignment_file])

    if tree_file is not None:
        cmd.extend(["--tree",tree_file])

    if model is not None:
        cmd.extend(["--model",model])

    # seed argument is overloaded. Interpret based on type
    if seed is not None:
        if type(seed) is bool:
            if seed:
                cmd.extend(["--seed",interface.gen_seed()])
        elif type(seed) is int:
            cmd.extend(["--seed",f"{seed:d}"])
        elif type(seed) is str:

            try:
                int(seed)
            except ValueError:
                err = f"seed {seed} could not be interpreted as an int\n"
                raise ValueError(err)

            cmd.extend(["--seed",seed])
        else:
            err = "seed must be True/False, int, or string representation of int\n"
            raise ValueError(err)

    num_threads = threads.get_num_threads(num_threads)
    cmd.extend(["--threads",f"{num_threads:d}"])

    # Put on any custom args
    if other_args is not None:
        for a in other_args:
            cmd.append(a)

    # If logging to standard out, get log file name
    log_file = None
    if log_to_stdout:
        log_file = "alignment.phy.raxml.log"

    # Run job
    interface.launch(cmd,dir_name,log_file)

    return " ".join(cmd)
