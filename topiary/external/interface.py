__description__ = \
"""
Generic interface for wrapping external programs.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-05-17"

import topiary

import pandas as pd

import subprocess, os, sys, time, random, string, shutil, copy, json
import multiprocessing as mp

def gen_seed():
    """
    Generate a random string of 10 integers and return as a string.
    """

    return "".join([f"{random.choice(range(10)):d}" for _ in range(10)])

def create_new_dir(dir_name=None,overwrite=False):
    """
    Create a new directory.

    dir_name: if specified, name the directory this
    overwrite: overwrite existing directory

    returns name of created directory
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

    return dir_name

def copy_input_file(input_file,
                    dir_name,
                    file_name=None,
                    make_input_dir=True,
                    df=None):
    """
    copy an input file into a directory in a stereotyped way.

    If make_input_dir is specified, copy input_file into dir_name/input,
    creating input if necessary.  If make_input_dir is not specified,
    copy in the file as dir_name/{input_file}.

    input_file: file to copy in
    dir_name: copy into dir_name
    file_name: what to call file in new directory. If None, use same name.
    make_input_dir: (bool) make input directory input or not.
    df: topiary data frame

    returns name of copied file
    """

    if file_name is None:
        file_name = os.path.split(input_file)[-1]
    file_alone = os.path.split(file_name)[-1]

    # If we are putting this into an input subdirectory...
    if make_input_dir:
        input_dir = os.path.join(dir_name,"input")
        if not os.path.exists(input_dir):
            os.mkdir(input_dir)
        file_alone = os.path.join("input",file_alone)

    # Copy in file
    out_file = os.path.join(dir_name,file_alone)
    shutil.copy(input_file,out_file)

    return file_alone


def prep_calc(previous_dir=None,
              df=None,
              model=None,
              tree_file=None,
              other_files=[],
              output=None,
              overwrite=False,
              output_base="calculation"):
    """
    Prepare a calculation, organizing input files and creating an output
    directory.

    previous_dir: directory containing previous calculation. prep_calc will
                  grab the the csv, model, and tree from the previous run.
    These arguments will override anything pulled out from the previous_dir.
    df: topiary data frame or .csv file written out from a topiary data frame
    model: phylogenetic model as string (e.g. JTT, LG+G8)
    tree_file: newick tree file
    other_files: other files besides the df and tree that are needed
    output: output directory. If not specified, create an output directory with
            form "{output_base}_randomletters"
    overwrite: whether or not to overwrite existing (default is False)
    output_base: base to assign the directory if no output is specified.
    """

    # -------------------------------------------------------------------------
    # Load in information from the previous calculation
    previous = {}
    if previous_dir is not None:
        previous = read_previous_run_dir(previous_dir)

    # -------------------------------------------------------------------------
    # Create output directory

    if output is None:
        rand = "".join([random.choice(string.ascii_letters) for _ in range(10)])
        output = f"{output_base}_{rand}"

    dir_name = create_new_dir(dir_name=output,overwrite=overwrite)

    # -------------------------------------------------------------------------
    # df (dataframe)

    # If passed df is None, try to grab previous.
    if df is None:
        try:
            df = previous["df"]
        except KeyError:
            pass

    # Read dataframe
    csv_file = None
    if df is not None:

        # Read dataframe
        if type(df) is pd.DataFrame:
            df = topiary.util.check_topiary_dataframe(df)

        # Read csv
        elif type(df) is str:
            csv_file = copy.deepcopy(df)
            df = topiary.read_dataframe(csv_file)
            csv_file = copy_input_file(csv_file,dir_name,make_input_dir=True)

        else:
            err = "df must be a pandas data frame or csv file that can be read as\n"
            err += "a pandas data frame\n"
            raise ValueError(err)

    # -------------------------------------------------------------------------
    # Tree file

    # If passed tree_file is None, try to grab previous
    if tree_file is None:
        try:
            tree_file = previous["tree_file"]
        except KeyError:
            pass

    # Copy in tree file
    if tree_file is not None:
        tree_file = copy_input_file(tree_file,dir_name,make_input_dir=True)

    # -------------------------------------------------------------------------
    # Model

    # If passed model is None, try to grab previous
    if model is None:
        try:
            model = previous["model"]
        except KeyError:
            pass

    if model is not None:
        if type(model) is not str:
            err = f"\nmodel '{model}' not recognized. Should be a string.\n"
            raise ValueError(err)

    # -------------------------------------------------------------------------
    # Copy remaining files into input directory

    final_files = []
    for f in other_files:
        if f is None:
            final_files.append(None)
        else:
            final_files.append(copy_input_file(f,
                                               dir_name,
                                               make_input_dir=True,
                                               df=df))

    # Move into the output directory
    starting_dir = os.getcwd()
    os.chdir(dir_name)

    # Make sure we actually generated the input directory at some point in the
    # steps above
    if not os.path.exists("input"):
        os.mkdir("input")

    # Write the alignment to the input directory
    alignment_file = os.path.join("input","alignment.phy")
    topiary.write_phy(df,out_file=alignment_file,seq_column="alignment",
                      write_only_keepers=True,clean_sequence=True)

    # write model to input directory
    if model is not None:
        f = open(os.path.join("input","model.txt"),"w")
        f.write(f"{model}\n")
        f.close()

    df_file = os.path.join("input","dataframe.csv")
    topiary.write_dataframe(df,out_file=df_file)

    out = {"df":df,
           "csv_file":csv_file,
           "tree_file":tree_file,
           "model":model,
           "alignment_file":alignment_file,
           "previous_dir":previous_dir,
           "starting_dir":starting_dir,
           "output":output,
           "other_files":final_files}

    return out

def _subproc_wrapper(cmd,stdout,queue):
    """
    Wrap the subprocess.run call to allow multithreading.

    args: args to pass to subprocess.run
    kwargs: kwargs to pass to subprocess.run
    queue: multiprocessing queue to catch return value
    """

    ret = subprocess.run(cmd,stdout=stdout)
    queue.put(ret)

def _follow_log_generator(f,queue):
    """
    Generator function that follows some file object (f) until a
    multiprocessing Queue (queue) is not empty. This is useful for following a
    log file being spit out by an external program whose return will be put
    into the queue.

    f: open file object
    queue: multiprocessing.Queue object
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

def launch(cmd,run_directory,log_file=None):
    """
    Launch an external command in a specific directory. Runs command on its own
    thread, allowing python to capture output from a file as standard output.

    cmd: subprocess style command list (i.e. ["ls","-al"])
    run_directory: directory in which to run the command. changes into that
                   directory, then returns out when process completes
    log_file: log file where command output will be stored. if specified, the
              output of the log file is captured and written to standard output
              (equivalent to `tail -f log_file`).

    returns None. Throws RuntimeError if the command terminates unexpectedly.
    """

    # Go into working directory
    cwd = os.getcwd()
    os.chdir(run_directory)

    # Print command
    full_cmd = " ".join(cmd)
    print(f"Running '{full_cmd}'",flush=True)

    # Launch as a multiprocessing process that will return its output to a
    # multiprocessing queue.
    queue = mp.Queue()
    main_process = mp.Process(target=_subproc_wrapper,
                              args=(cmd,subprocess.PIPE,queue))
    main_process.start()

    # If following a log
    if log_file is not None:

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

        err += "".join([line for line in ret.stdout.decode()])

        raise RuntimeError(err)

    # Leave working directory
    os.chdir(cwd)

def write_run_information(outdir,df,calc_type,model,cmd):
    """
    Write information from the run in a standard way.

    Parameters
    ----------
        outdir: output directory (i.e. ml-run/output)
        df: dataframe to write to dataframe.csv
        calc_type: calculation type (string)
        model: model (string)
        cmd: command invoked to raxml or generax

    Return
    ------
        None
    """

    # Dataframe
    topiary.write_dataframe(df,os.path.join(outdir,"dataframe.csv"))

    # Run parameters
    out_dict = {"calc_type":calc_type,
                "model":model,
                "cmd":cmd,
                "version":topiary.__version__}
    f = open(os.path.join(outdir,"run_parameters.json"),"w")
    json.dump(out_dict,f)
    f.close()


def read_previous_run_dir(previous_dir):
    """
    Load the df, tree, and run_parameters from a previous run directory.

    Parameters
    ----------
        previous_dir: directory of previous run (top level; meaning, ml-tree
                      not ml-tree/output)

    Return
    ------
        dictionary with df, tree_file, and parameters
    """

    previous = {}

    # Make sure previous_dir is a string.
    if type(previous_dir) is not str:
        err = f"\nprevious_dir '{previous_dir}' not recognized. Should be a string.\n"
        raise ValueError(err)

    # Look for output directory
    out_dir = os.path.join(previous_dir,"output")
    if not os.path.isdir(out_dir):
        err = f"\nCould not read previous directory '{previous_dir}'. This\n"
        err += "should be a previous topiary run directory that contains an\n"
        err += "'output' directory.\n"
        raise ValueError(err)

    # Grab the run parameters
    try:
        f = open(os.path.abspath(os.path.join(out_dir,"run_parameters.json")),'r')
        run_parameters = json.load(f)
        f.close()
    except FileNotFoundError:
        err = f"\nCould not read previous directory '{previous_dir}'. This\n"
        err += "directory should contain a file output/run_parameters.json.\n"
        err += "This file was not found.\n\n"
        raise ValueError(err)

    # copy in previous run parameters
    for k in run_parameters:
        previous[k] = run_parameters[k]

    # Try to grab dataframe
    df_file = os.path.abspath(os.path.join(out_dir,"dataframe.csv"))
    if os.path.exists(df_file):
        previous["df"] = topiary.read_dataframe(df_file)
    else:
        err = f"\nCould not read previous directory '{previous_dir}'. Does not\n"
        err += "contain the file output/dataframe.csv"
        raise ValueError(err)

    # Try to grab the tree file
    tree_file = os.path.abspath(os.path.join(out_dir,"tree.newick"))
    if os.path.exists(tree_file):
        previous["tree_file"] = tree_file

    return previous
