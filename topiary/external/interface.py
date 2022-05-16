
import topiary

import pandas as pd

import subprocess, os, sys, time, random, string, shutil, copy
import multiprocessing as mp

def gen_seed():
    """
    Generate a random string of 10 integers and return as a string.
    """

    return "".join([f"{random.choice(range(10)):d}" for _ in range(10)])

def create_new_dir(dir_name=None):
    """
    Create a new directory.

    dir_name: if specified, name the directory this

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

    # If we are putting this into an input subdirectory
    if make_input_dir:
        input_dir = os.path.join(dir_name,"input")
        if not os.path.exists(input_dir):
            os.mkdir(input_dir)
        file_alone = os.path.join("input",file_alone)

    # Copy in file, potentially converting from to uid
    out_file = os.path.join(dir_name,file_alone)
    shutil.copy(input_file,out_file)

    return file_alone

def _load_previous_dir(previous_dir):
    """
    Load the df, tree, and model from a previous run directory.

    previous_dir: output directory from previous run

    returns: dictionary with df, tree_file, and model keys (if each element found)
    """

    previous = {}

    # Make sure previous_dir is a string.
    if type(previous_dir) is not str:
        err = f"\nprevious_dir '{previous_dir}' not recognized. Should be a string.\n"
        raise ValueError(err)

    out_dir = os.path.join(previous_dir,"output")
    if not os.path.isdir(out_dir):
        err = f"\nCould not read previous directory '{previous_dir}'. This\n"
        err += "should be a previous topiary run directory that contains an\n"
        err += "'output' directory.\n"
        raise ValueError(err)

    # Try to grab dataframe
    df_file = os.path.abspath(os.path.join(out_dir,"dataframe.csv"))
    if os.path.exists(df_file):
        previous["df"] = topiary.read_dataframe(df_file)

    # Try to grab the tree file
    tree_file = os.path.abspath(os.path.join(out_dir,"tree.newick"))
    if os.path.exists(tree_file):
        previous["tree_file"] = tree_file

    # Try to grab the model
    model_file = os.path.abspath(os.path.join(out_dir,"model.txt"))
    if os.path.exists(model_file):

        f = open(model_file,'r')
        model = f.read().strip()
        f.close()

        previous["model"] = model

    return previous


def prep_calc(previous_dir=None,
              output=None,
              df=None,
              model=None,
              tree_file=None,
              other_files=[],
              output_base="calculation"):
    """
    Prepare a calculation, organizing input files and creating an output
    directory.

    previous_dir: directory containing previous calculation. prep_calc will
                  grab the the csv, model, and tree from the previous run.
    output: output directory. If not specified, create an output directory with
            form "{output_base}_randomletters"

    These arguments will override anything pulled out from the previous_dir.
    df: topiary data frame or .csv file written out from a topiary data frame
    model: phylogenetic model as string (e.g. JTT, LG+G8)
    tree_file: newick tree file
    other_files: other files besides the df and tree that are needed

    output_base: base to assign the directory if no output is specified.
    """

    # -------------------------------------------------------------------------
    # Load in information from the previous calculation
    previous = {}
    if previous_dir is not None:
        previous = _load_previous_dir(previous_dir)

    # -------------------------------------------------------------------------
    # Create output directory

    if output is None:
        rand = "".join([random.choice(string.ascii_letters) for _ in range(10)])
        output = f"{output_base}_{rand}"

    dir_name = create_new_dir(dir_name=output)

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
    alignment_file = os.path.join("input","alignment")
    topiary.write_phy(df,out_file=alignment_file,seq_column="alignment",
                      write_only_keepers=True,clean_sequence=True)

    out = {"df":df,
           "csv_file":csv_file,
           "tree_file":tree_file,
           "model":model,
           "alignment_file":alignment_file,
           "previous_dir":previous_dir,
           "starting_dir":starting_dir,
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

def _follow_log_generator(f,p):
    """
    Generator function that follows some file object (f) until some
    multiprocessing Process (p) is not longer alive. This is useful for
    following a log file being spit out by an external program running on p.

    f: open file object
    p: multiprocessing.Process object
    """

    # start infinite loop
    while p.is_alive():
        # read last line of file
        line = f.readline()
        # sleep if file hasn't been updated
        if not line:
            time.sleep(0.1)
            continue

        yield line

def launch(cmd,run_directory,log_file=None):

    # Go into working directory
    cwd = os.getcwd()
    os.chdir(run_directory)

    # Print command
    full_cmd = " ".join(cmd)
    print(f"Running '{full_cmd}'",flush=True)

    # Launch generax as a multiprocessing process dumping its output to a
    # multiprocessing queue.
    queue = mp.Queue()
    main_process = mp.Process(target=_subproc_wrapper,
                              args=(cmd,subprocess.PIPE,queue))
    main_process.start()

    # If dumping log
    if log_file is not None:

        # While main process is running
        while main_process.is_alive():

            # If queue is empty, raxml job hasn't finished yet
            if not queue.empty():
                break

            # Try to open log every second
            try:
                f = open(log_file,"r")
                # Use follow generator function to catch lines as the come out
                for line in _follow_log_generator(f,main_process):
                    sys.stdout.write(line)
                    sys.stdout.flush()
                f.close()
            except FileNotFoundError:
                time.sleep(1)

    # Wait for main process to complete and get return
    main_process.join()
    ret = queue.get()

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
