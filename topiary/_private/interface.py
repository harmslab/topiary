"""
Generic interface for wrapping external programs.
"""

import topiary
from topiary._private import check

import pandas as pd

import subprocess
import os
import sys
import time
import random
import string
import shutil
import copy
import json
import glob
import multiprocessing as mp

class MockTqdm():
    """
    Fake tqdm progress bar so we don't have to show a status bar if we don't
    want to. Can be substituted wherever we would use tqdm (i.e.
    tqdm(range(10)) --> MockTqdm(range(10)).
    """

    def __init__(self,*args,**kwargs):
        pass
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        pass

    def update(self,*args,**kwargs):
        pass


def gen_seed():
    """
    Generate a random string of 10 integers and return as a string.

    Return
    ------
    seed : str
        10-digit integer as a string
    """

    return "".join([f"{random.choice(range(10)):d}" for _ in range(10)])

def create_new_dir(dir_name=None,overwrite=False):
    """
    Create a new directory.

    Parameters
    ----------
    dir_name : str
        if specified, name the directory this
    overwrite : bool, default=False
        overwrite existing directory

    Return
    ------
    directory : str
        name of created directory
    """

    # if dir_name is not specified, build it in a stereotyped fashion
    if dir_name is None:
        rand_name = "".join([random.choice(string.ascii_letters)
                              for _ in range(10)])
        dir_base = "calculation"

        dir_name = f"{dir_base}_{rand_name}"

    # If directory already exists, throw error
    if os.path.exists(dir_name):
        if overwrite:
            shutil.rmtree(dir_name)
        else:
            err = f"{dir_name} already exists.\n"
            raise FileExistsError(err)

    # Make directory
    os.mkdir(dir_name)

    return os.path.abspath(dir_name)

def copy_input_file(input_file,
                    dir_name,
                    file_name=None):
    """
    Copy an input file into a directory.

    Parameters
    ----------
    input_file : str
        file to copy in
    dir_name : str
        copy into dir_name
    file_name : str, optional
        what to call file in new directory. If None, use same name.

    Return
    ------
    file : str
        name of copied file
    """

    if not os.path.isfile(input_file):
        err = f"input_file {input_file} does not exist\n"
        raise FileNotFoundError(err)

    if not os.path.isdir(dir_name):
        err = f"dir_name {dir_name} is not a directory\n"
        raise FileNotFoundError(err)

    if file_name is None:
        file_name = input_file

    file_alone = os.path.basename(file_name)

    # Copy in file
    out_file = os.path.join(dir_name,file_alone)
    shutil.copy(input_file,out_file)

    return os.path.abspath(out_file)

def _follow_log_subproc_wrapper(cmd,stdout,queue):
    """
    Wrap the subprocess.run call to allow multithreaded following of a log file.

    Parameters
    ----------
    cmd : list
        args to pass to subprocess.run
    stdout : subprocess.PIPE
        where to write standard output
    queue : multiprocessing.Queue
        multiprocessing queue to catch return value

    Return
    ------
    None
    """

    ret = subprocess.run(cmd,stdout=stdout)
    queue.put(ret)

def _follow_log_generator(f,queue):
    """
    Generator function that follows some file object (f) until a
    multiprocessing Queue (queue) is not empty. This is useful for following a
    log file being spit out by an external program whose return will be put
    into the queue.

    Parameters
    ----------
    f : open file object
        log file to read
    queue: multiprocessing.Queue object
        multiprocessing queue to append output from file

    Return
    ------
    None
    """

    # start infinite loop
    found_line = False
    while queue.empty():

        # read last line of file
        line = f.readline()

        # sleep if file hasn't been updated
        if not line:
            time.sleep(0.1)
            continue

        found_line = True

        yield line

    # If we found a line at one point and the job has finished, keep reading
    # lines until we run out -- that is, get trailing data written to log file
    # after job is done
    counter = 0
    if found_line and not queue.empty():
        while True:
            counter += 1
            line = f.readline()
            if line is None or line == "":
                return

            yield line

            if counter > 200:
                break

def launch(cmd,
           run_directory,
           log_file=None,
           suppress_output=False,
           write_to_script=None):
    """
    Launch an external command in a specific directory. If log_file is
    specified, runs command on its own thread, allowing python to capture output
    from a file as standard output.

    Parameters
    ----------
    cmd : list
        subprocess style command list (i.e. ["ls","-al"])
    run_directory : str
        directory in which to run the command. changes into that directory, then
        returns back from that directory when process completes
    log_file : str, optional
        log file where command output will be stored. if specified, the
        output of the log file is captured and written to standard output
        (equivalent to `tail -f log_file`).
    suppress_output : bool, default=False
        whether or not to capture (and not return) stdout and stderr. Ignored
        if log_file is specified.
    write_to_script : str, optional
        instead of running the command, write out the command to the script file
        in the run directory. this can then be invoked later by something like
        :code:`bash script_file`.

    Returns
    -------
    None

    Raises
    ------
    RuntimeError :
        If the command terminates unexpectedly.
    """

    # Go into working directory
    cwd = os.getcwd()
    os.chdir(run_directory)

    # Print command
    full_cmd = " ".join(cmd)
    if not suppress_output:
        print(f"Running '{full_cmd}'",flush=True)

    # If write_to_Script is specified, write that out instead of actually
    # running. Return to starting directory and return None.
    if write_to_script is not None:

        f = open(str(write_to_script),"w")
        f.write(full_cmd)
        if suppress_output:
            f.write(" &> topiary.log")
        f.write("\n")
        f.close()

        os.chdir(cwd)
        return None

    # If no log file specified, run directly
    if log_file is None:
        if suppress_output:
            capture_output = True
        else:
            capture_output = False
        ret = subprocess.run(cmd,capture_output=capture_output)

    # Otherwise, run on it's own thread and capture output to standard out
    else:

        # Launch as a multiprocessing process that will return its output to a
        # multiprocessing queue.
        queue = mp.Queue()
        main_process = mp.Process(target=_follow_log_subproc_wrapper,
                                  args=(cmd,subprocess.PIPE,queue))
        main_process.start()

        # If queue is not empty, the job has finished and put its return value
        # into the queue
        while queue.empty():

            # Try to open log every second
            try:
                f = open(log_file,"r")

                # Use follow generator function to catch lines as the come out
                for line in _follow_log_generator(f,queue):
                    print(line,flush=True,end="")
                f.close()

            except FileNotFoundError:
                time.sleep(1)

        # Wait for main process to complete and get return
        ret = queue.get()
        main_process.join()

    # Check for error on return
    if ret.returncode != 0:
        err = f"ERROR: {cmd[0]} returned {ret.returncode}\n\n"
        err += "------------------------------------------------------------\n"
        err += f" {cmd[0]} output \n"
        err += "------------------------------------------------------------\n"
        err += "\n\n"

        if ret.stdout is not None:
            err += "".join([line for line in ret.stdout.decode()])

        raise RuntimeError(err)

    # Leave working directory
    os.chdir(cwd)
