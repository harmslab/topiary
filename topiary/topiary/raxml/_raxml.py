__description__ = \
"""
Core functions for wrapping raxml-ng. 
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2021-07-22"

# raxml binary to use it not specified by user
RAXML_BINARY = "raxml-ng.dev"

import topiary

import pandas as pd

import subprocess, os, sys, time, random, string, shutil
import multiprocessing as mp


def gen_seed():
    """
    Generate a random string of 10 integers and return as a string for passing
    to raxml.
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
        dir_base = os.path.split(RAXML_BINARY)[-1]

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
                     pretty_to_uid=False,
                     df=None):
    """
    copy an input file into a directory in a stereotyped way.

    If make_input_dir is specified, copy input_file into dir_name/00_input,
    creating 00_input if necessary.  If make_input_dir is not specified,
    copy in the file as dir_name/{input_file}.

    input_file: file to copy in
    dir_name: copy into dir_name
    file_name: what to call file in new directory. If none, use same name.
    make_input_dir: (bool) make input directory 00_input or not.
    pretty_to_uid: (bool) convert file to uid
    df: topiary data frame (required if pretty_to_uid is True)

    returns name of copied file
    """


    if file_name is None:
        file_name = os.path.split(input_file)[-1]
    file_alone = os.path.split(file_name)[-1]

    # If we are putting this into an input subdirectory
    if make_input_dir:
        input_dir = os.path.join(dir_name,"00_input")
        if not os.path.exists(input_dir):
            os.mkdir(input_dir)
        file_alone = os.path.join("00_input",file_alone)

    # Copy in file, potentially converting from to uid
    out_file = os.path.join(dir_name,file_alone)
    if pretty_to_uid:
        if df is None:
            err = "you must specify df if you set pretty_to_uid = True\n"
            raise ValueError(err)
        topiary.util.pretty_to_uid(df,input_file,out_file=out_file)
    else:
        shutil.copy(input_file,out_file)

    return file_alone

def prep_calc(df=None,
              output=None,
              other_files=[],
              output_base="raxml-calc"):
    """
    Prepare a raxml calculation, generating raxml-readable input files,
    creating an output directory, and moving into the output directory.

    df: topiary data frame or .csv file written out from a topiary data frame
    output: output directory. If not specified, create an output directory with
            form "{output_base}_randomletters"
    other_file: other files besides the df needed (usually tree files)
    output_base: base to assign the directory if no output is specified.
    """

    # Create output directory
    if output is None:
        rand = "".join([random.choice(string.ascii_letters) for _ in range(10)])
        output = f"{output_base}_{rand}"

    dir_name = _create_new_dir(dir_name=output)

    # Deal with data frame
    if df is not None:

        # Read/parse csv file.
        csv_file = None
        if type(df) is pd.DataFrame:
            pass
        elif type(df) is str:
            csv_file = df
            df = pd.read_csv(csv_file)
        else:
            err = "df must be a pandas data frame or csv file that can be read as\n"
            err += "a pandas data frame\n"
            raise ValueError(err)
    else:
        df = None
        csv_file = None

    if csv_file is not None:
        csv_file = _copy_input_file(csv_file,dir_name,make_input_dir=True)

    # Copy files into input directory. This will put them in 00_input and keep
    # their original filenames so we have some notion of where they came from.
    final_files = []
    for f in other_files:
        if f is None:
            final_files.append(None)
        else:
            final_files.append(_copy_input_file(f,
                                                dir_name,
                                                make_input_dir=True,
                                                pretty_to_uid=True,
                                                df=df))

    # Move into the output directory
    starting_dir = os.getcwd()
    os.chdir(dir_name)

    # Write out alignment
    if not os.path.exists("00_input"):
        os.mkdir("00_input")
    alignment_file = os.path.join("00_input","alignment")
    topiary.write_phy(df,out_file=alignment_file,seq_column="alignment",
                      write_only_keepers=True,clean_sequence=True)

    out = {"df":df,
           "csv_file":csv_file,
           "alignment_file":alignment_file,
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
    dir_name = _create_new_dir(dir_name=dir_name)

    # Copy alignment and tree files into the directory (if specified)
    if alignment_file is not None:
        alignment_file = _copy_input_file(alignment_file,
                                          dir_name,
                                          file_name="alignment",
                                          make_input_dir=False)
    if tree_file is not None:
        tree_file = _copy_input_file(tree_file,
                                     dir_name,
                                     file_name="tree",
                                     make_input_dir=False)

    # Go into working directory
    cwd = os.getcwd()
    os.chdir(dir_name)

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
            cmd.extend(["--seed",_gen_seed()])
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

    # Construct command and dump to std out
    full_cmd = " ".join(cmd)
    print(f"Running '{full_cmd}'")
    sys.stdout.flush()

    # Launch raxml as a multiprocessing process dumpint its output to a
    # multiprocessing queue.
    queue = mp.Queue()
    main_process = mp.Process(target=_subproc_wrapper,
                              args=(cmd,subprocess.PIPE,queue))
    main_process.start()

    # If dumping log
    if log_to_stdout:

        # While main process is running
        while main_process.is_alive():

            # Try to open log every second
            try:
                f = open("alignment.raxml.log","r")
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
        err = f"ERROR: raxml returned {ret.returncode}\n\n"
        err += "------------------------------------------------------------\n"
        err += " raxml output \n"
        err += "------------------------------------------------------------\n"
        err += "\n\n"

        err += "".join([line for line in ret.stdout.decode()])

        raise RuntimeError(err)

    # Leave working directory
    os.chdir(cwd)
