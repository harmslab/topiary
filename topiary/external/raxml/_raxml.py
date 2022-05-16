__description__ = \
"""
Core functions for wrapping raxml-ng.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2021-07-22"

# raxml binary to use it not specified by user
RAXML_BINARY = "raxml-ng"

import topiary
from topiary.external import interface

import pandas as pd

def run_raxml(algorithm=None,
              alignment_file=None,
              tree_file=None,
              model=None,
              dir_name=None,
              seed=None,
              threads=1,
              raxml_binary=RAXML_BINARY,
              log_to_stdout=True,
              other_args=[]):
    """
    Run raxml. Creates a working directory, copies in the relevant files, runs
    there, and then returns to the previous directory.

    algorithm: algorithm to run (--all, --ancestral, etc.)
    alignment_file: alignment file in .phy format (passed via --msa)
    tree_file: tree file in .newick format (passed via --tree)
    model: model in format recognized by --model
    dir_name: If specified, this will be the name of the working directory.
    seed: true/false, int, or str. If true, pass a randomly generated seed to
          raxml. If int or str, use that as the seed. (passed via --seed)
    threads: number of threads to use (passed via --threads)
    raxml_binary: raxml binary to use
    log_to_stdout: capture log and write to std out.
    other_args: list of arguments to pass to raxml
    """

    # Create directory in which to do calculation
    dir_name = interface.create_new_dir(dir_name=dir_name)

    # Copy alignment and tree files into the directory (if specified)
    if alignment_file is not None:
        alignment_file = interface.copy_input_file(alignment_file,
                                                   dir_name,
                                                   file_name="alignment",
                                                   make_input_dir=False)
    if tree_file is not None:
        tree_file = interface.copy_input_file(tree_file,
                                              dir_name,
                                              file_name="tree",
                                              make_input_dir=False)

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

    cmd.extend(["--threads",f"{threads:d}"])

    # Put on any custom args
    for a in other_args:
        cmd.append(a)

    # If logging to standard out, get log file name
    log_file = None
    if log_to_stdout:
        log_file = "alignment.raxml.log"

    # Run job
    interface.launch(cmd,dir_name,log_file)
