__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
BLAST against a local database.
"""

from . base import read_blast_xml

import Bio.Blast.Applications as apps

import pandas as pd
from tqdm.auto import tqdm

import sys, os, string, random, subprocess, copy
import warnings
import multiprocessing as mp

def _blast_thread(args):
    """
    Run reverse blast on a thread. Should only be called via _blast_thread_manager.

    takes args which are interpreted as (sequence_list,i,blast_kwargs,
                                         blast_function,keep_tmp,queue)

    sequence_list: list of sequences as strings
    i: sequence to grab
    blast_kwargs: keyword arguments to pass to blast call
    blast_function: blast function to call
    keep_tmp: whether or not to keep temporary files
    queue: multiprocessing queue for storing results

    puts resulting hit as a pandas dataframe into the queue
    """

    # parse args
    sequence_list = args[0]
    i = args[1]
    blast_function = args[2]
    blast_kwargs = args[3]
    keep_tmp = args[4]
    queue = args[5]

    # make a 10-character random string for temporary files
    tmp_file_root = "".join([random.choice(string.ascii_letters) for i in range(10)])
    input_file = "{}_blast_in.fasta".format(tmp_file_root)
    out_file = "{}_blast_out.xml".format(tmp_file_root)

    f = open(input_file,'w')
    f.write("".join(f">sequence{i}\n{sequence_list[i]}"))
    f.close()

    blast_kwargs = copy.deepcopy(blast_kwargs)
    blast_kwargs["query"] = input_file
    blast_kwargs["out"] = out_file

    blast_function(**blast_kwargs)()

    # Parse output
    try:
        out_df = read_blast_xml(out_file)
    except FileNotFoundError:
        err = "\nLocal blast failed on sequence:\n"
        err += f"    '{sequence_list[i]}'\n\n"
        raise RuntimeError(err)

    # Delete temporary files
    if not keep_tmp:
        os.remove(input_file)
        os.remove(out_file)

    # update queue
    queue.put((i,out_df))


def _blast_thread_manager(sequence_list,
                          blast_function,
                          blast_kwargs,
                          keep_tmp=False,
                          num_threads=-1):
    """
    Run a bunch of blast jobs in a mulithreaded fashion. Should only be called
    by local_blast.

    sequence_list: list of sequences as strings
    blast_function: blast function to use
    blast_kwargs: keyword arguments to pass to blast call
    keep_tmp: whether or not to keep temporary files
    num_threads: number of threads to use. if -1, use all available.

    Returns a list of dataframes with blast results ordered by the input
    sequence list.
    """

    print("Performing blast...")
    sys.stdout.flush()

    # Figure out number of threads to use
    if num_threads < 0:
        try:
            num_threads = mp.cpu_count()
        except NotImplementedError:
            num_threads = os.cpu_count()
            if num_threads is None:
                warnings.warn("Could not determine number of cpus. Using single thread.\n")
                num_threads = 1

    # queue will hold results from each run.
    queue = mp.Manager().Queue()
    with mp.Pool(num_threads) as pool:

        # This is a bit obscure. Build a list of args to pass to the pool. Each
        # tuple of args matches the args in _reverse_blast_thread.
        # all_args has all len(df) reverse blast runs we want to do.
        all_args = []
        for i in range(len(sequence_list)):
            all_args.append((sequence_list,
                             i,
                             blast_function,
                             blast_kwargs,
                             keep_tmp,
                             queue))

        # Black magic. pool.imap() runs a function on elements in iterable,
        # filling threads as each job finishes. (Calls _blast_thread
        # on every args tuple in all_args). tqdm gives us a status bar.
        # By wrapping pool.imap iterator in tqdm, we get a status bar that
        # updates as each thread finishes.
        list(tqdm(pool.imap(_blast_thread,all_args),total=len(all_args)))
        
    # Get results out of the queue.
    results = []
    while not queue.empty():
        results.append(queue.get())

    # Sort results
    results.sort()
    hits = [r[1] for r in results]

    print("Success.")
    sys.stdout.flush()

    return hits

def local_blast(sequence,
                db,
                blast_program="blastp",
                hitlist_size=100,
                e_value_cutoff=0.01,
                gapcosts=(11,1),
                keep_tmp=False,
                num_threads=-1,
                **kwargs):
    """
    Perform a blast query against a local blast data base.

    ssequence: sequence as a string OR list of string sequences
    db: name of local blast database
    blast_program: NCBI blast program to use (blastp, tblastn, etc.)
    hitlist_size: download only the top hitlist_size hits
    e_value_cutoff: only take hits with e_value better than e_value_cutoff
    gapcosts: BLAST gapcosts (length 2 tuple of ints)
    keep_tmp: whether or not to keep temporary blast output
    num_threads: number of threads to use. if -1, use all available.
    kwargs: extra keyword arguments are passed directly to apps.NcbiblastXXXCommandline,
            overriding anything constructed above.
    """

    recognized_functions = {"blastp":apps.NcbiblastpCommandline,
                            "blastn":apps.NcbiblastnCommandline,
                            "blastx":apps.NcbiblastxCommandline,
                            "tblastn":apps.NcbitblastnCommandline,
                            "tblastx":apps.NcbitblastxCommandline,
                            "psiblast":apps.NcbipsiblastCommandline,
                            "rpsblast":apps.NcbirpsblastCommandline,
                            "rpstblastn":apps.NcbirpstblastnCommandline,
                            "deltablast":apps.NcbideltablastCommandline}

    if not os.path.exists(f"{db}.psq"):
        err = f"{db} not found!\n"
        raise RuntimeError(err)

    # Make sure we can actually run the local blasting
    try:
        blast_function = recognized_functions[blast_program]
    except KeyError:
        err = "\nblast_program '{}' not recognized.\n\nAllowed programs:\n".format(blast_program)
        for k in recognized_functions.keys():
            err += "    {}\n".format(k)
        raise ValueError(err)

    # Make sure the blast program is installed.
    try:
        ret = subprocess.run([blast_program],stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        err = f"\nBLAST program {_local_blast} is not in path. Is it installed?\n\n"
        raise FileNotFoundError(err)

    # Parse the sequence input. Should either be a list of strings or a single
    # string sequence. Create sequence_list which is a list of sequences (even
    # if only one)
    if hasattr(sequence,"__iter__"):

        # If sequence is not a string, make a single fasta-formatted string
        # from contents
        if type(sequence) is not str:

            # pandas dataframes will pass below because enumeration goes over
            # columns, which are strings. Manually check for dataframe...
            if type(sequence) is pd.DataFrame:
                err = "\nsequence must either be a string or list of strings, one for each sequence\n\n"
                raise ValueError(err)

            sequence_list = []
            for i, s in enumerate(sequence):
                if type(s) is not str:
                    err = "\nsequence must either be a string or list of strings, one for each sequence\n\n"
                    raise ValueError(err)
                else:
                    sequence_list.append(s)
        else:
            sequence_list = [sequence]
    else:
        err = "\nsequence must either be a string or list of strings, one for each sequence\n\n"
        raise ValueError(err)

    # Parse gap penalties
    try:
        bad = True
        if len(gapcosts) == 2:
            if type(gapcosts[0]) is int and type(gapcosts[1]) is int:
                bad = False
        if bad:
            raise ValueError
    except (TypeError,ValueError):
        err = "\ngapcosts should be a length-two tuple of ints\n\n"
        raise ValueError(err)

    gaps = '{} {}'.format(*gapcosts)

    # Construct keyword arguments to pass to function
    blast_kwargs = {"cmd":blast_program,
                    "db":db,
                    "outfmt":5,
                    "max_target_seqs":hitlist_size,
                    "threshold":e_value_cutoff,
                    "gapopen":gapcosts[0],
                    "gapextend":gapcosts[1]}
    for k in kwargs:
        blast_kwargs[k] = kwargs[k]

    # Run multi-threaded local blast
    out_df = _blast_thread_manager(sequence_list,
                                   blast_function,
                                   blast_kwargs,
                                   keep_tmp,
                                   num_threads)

    # If there was only one blast sequence, return it instead of a list
    if len(out_df) == 1:
        out_df = out_df[0]

    return out_df
