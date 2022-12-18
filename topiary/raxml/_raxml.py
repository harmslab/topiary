"""
Run raxml. Creates a working directory, copies in the relevant files, runs
there, and then returns to the previous directory.
"""

# raxml binary to use if not specified by user
RAXML_BINARY = "raxml-ng"

import topiary
from topiary._private import interface
from topiary._private import threads
from topiary._private import Supervisor
from topiary._private import check

import os
import shutil

def run_raxml(run_directory,
              algorithm=None,
              alignment_file=None,
              tree_file=None,
              model=None,
              seed=None,
              log_to_stdout=True,
              suppress_output=False,
              other_args=None,
              other_files=None,
              write_to_script=None,
              supervisor=None,
              num_threads=-1,
              raxml_binary=RAXML_BINARY):
    """
    Run raxml. Creates a directory, copies in the relevant files, runs there,
    and then returns to the previous directory.

    Parameters
    ----------
    run_directory : str
        Name of the working directory.
    algorithm : str
        algorithm to run (--all, --ancestral, etc.)
    alignment_file : str
        alignment file in .phy format (passed via --msa)
    tree_file : str
        tree file in .newick format (passed via --tree)
    model : str
        model in format recognized by --model
    seed : bool,int,str
        If true, pass a randomly generated seed to raxml. If int or str, use
        that as the seed. (passed via --seed)
    log_to_stdout : bool, default=True
        capture log and write to std out.
    suppress_output : bool, default=False
        suppress output entirely. (ignored if log_to_stdout is True)
    other_args : list-like, optional
        list of arguments to pass to raxml
    other_files : list-like, optional
        list of files to copy into working directory (besides tree_file and
        alignment_file)
    write_to_script : str, optional
        instead of running the command, write out the command to the script file
        in the run directory. this can then be invoked later by something like
        :code:`bash script_file`.
    supervisor : Supervisor, optional
        instance of Supervisor to record when we start the job
    num_threads : int, default=-1
        number of threads (passed via --threads). if -1, use all available.
    raxml_binary : str, default=RAXML_BINARY
        raxml binary to use

    Return
    ------
    raxml_command : string
        command passed to raxml-ng as a string
    """

    if write_to_script is not None:
        if not issubclass(type(write_to_script),str):
            if supervisor is not None:
                supervisor.finalize(successful=False)
            err = "write_to_script should be None or a string giving script name\n"
            raise ValueError(err)

    # If the run_directory already exists...
    if os.path.exists(run_directory):
        if supervisor is not None:
            supervisor.finalize(successful=False)
        err = f"run_directory '{run_directory}' already exists\n"
        raise FileExistsError(err)

    # Make run directory
    os.mkdir(run_directory)

    # Copy alignment and tree files into the directory (if specified)
    if alignment_file is not None:
        shutil.copy(alignment_file,os.path.join(run_directory,"alignment.phy"))

    if tree_file is not None:
        shutil.copy(tree_file,os.path.join(run_directory,"tree.newick"))

    # Copy in any other required files, if requested
    if other_files is not None:
        for i in range(len(other_files)):
            file_name = os.path.basename(other_files[i])
            shutil.copy(other_files[i],os.path.join(run_directory,file_name))

    # Build a command list. Put in full path to raxml_binary
    abs_path_raxml_binary = shutil.which(raxml_binary)
    if abs_path_raxml_binary is None:
        if supervisor is not None:
            supervisor.finalize(successful=False)
        err = f"\nraxml_binary '{raxml_binary}' could not be found in the PATH\n\n"
        raise FileNotFoundError(err)

    cmd = [abs_path_raxml_binary]

    if algorithm is not None:
        cmd.append(algorithm)

    if alignment_file is not None:
        cmd.extend(["--msa","alignment.phy"])

    if tree_file is not None:
        cmd.extend(["--tree","tree.newick"])

    if model is not None:
        cmd.extend(["--model",model])

    # seed argument is overloaded. Interpret based on type
    if seed is not None:

        # If bool and True, make the seed
        try:
            seed = check.check_bool(seed)
            if seed:
                seed = interface.gen_seed()
            else:
                seed = 0
        except ValueError:
            pass

        # Make sure the seed -- whether passed in or generated above -- is
        # actually an int.
        try:
            seed = check.check_int(seed,minimum_allowed=0)
        except ValueError:
            if supervisor is not None:
                supervisor.finalize(successful=False)
            err = f"seed '{seed}' invalid. must be True/False or int > 0\n"
            raise ValueError(err)

        # If we have a seed > 0, append to command
        if seed > 0:
            cmd.extend(["--seed",f"{seed:d}"])

    # Figure out how to treat threads
    try:
        num_threads = threads.get_num_threads(num_threads)
    except ValueError as e:
        if supervisor is not None:
            supervisor.finalize(successful=False)
        raise ValueError from e


    if algorithm in ["--all","--search"]:
        threads_arg = "auto{" + f"{num_threads:d}" + "}"
    else:
        threads_arg = f"{num_threads:d}"

    cmd.extend(["--threads",threads_arg])

    # Put on any custom args
    if other_args is not None:
        for a in other_args:
            cmd.append(a)

    # If logging to standard out, get log file name
    log_file = None
    if log_to_stdout:
        log_file = "alignment.phy.raxml.log"

    if supervisor is not None:
        supervisor.event("Launching raxml-ng",
                         cmd=cmd,
                         num_threads=num_threads)

    # Run job
    try:
        interface.launch(cmd,
                         run_directory=run_directory,
                         log_file=log_file,
                         suppress_output=suppress_output,
                         write_to_script=write_to_script)
    except interface.WrappedFunctionException as e:
        if supervisor is not None:
            supervisor.finalize(successful=False)
        raise RuntimeError from e

    return " ".join(cmd)
