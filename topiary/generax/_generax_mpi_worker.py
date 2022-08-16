"""
MPI-aware worker script for generax. This should be called by the
reoncile_bootstrap function.
"""

import topiary
from mpi4py import MPI
import time
import os
import sys
import subprocess
import pathlib
import glob

def main(argv=None):
    """
    Run a generax worker calculation, checking for and avoiding collisions with
    other workers. Also gives complete restart capability to reconcile bootstrap
    calculation. Generally invoked from a command line with mpirun.
    """

    # Get command line arguments
    if argv is None:
        argv = sys.argv[1:]

    # Parse command line arguments
    try:
        replicates_dir = argv[0]
        run_id = argv[1]
    except (IndexError,):
        err = "\nreplicates_dir and run_id must be specified on the command line\n"
        raise ValueError(err)

    # Get MPI information (rank, particularly)
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Change into replicates_directory and get sorted list of bootstrap
    # replicates.
    os.chdir(replicates_dir)
    dirs = [d for d in os.listdir(".") if os.path.isdir(d)]
    dirs.sort()

    # Filenames for workers staking claims or declaring completion will start
    # with ...
    claim_base = f"claimed_run-{run_id}"
    complete_base = f"completed"
    skip_base = f"skipped"

    # Go through each directory
    for d in dirs:

        # If this directory has already been run
        completed_files = glob.glob(os.path.join(d,f"{complete_base}*"))
        if len(completed_files) > 0:
            continue

        # See if another worker has claimed this directory
        claim_files  = glob.glob(os.path.join(d,f"{claim_base}*"))
        if len(claim_files) > 0:
            continue

        # See if another worker has claimed this directory
        skip_files = glob.glob(os.path.join(d,f"{skip_base}*"))
        if len(skip_files) > 0:
            continue

        # Stake a claim
        os.chdir(d)
        my_claim_file = f"{claim_base}_by-{rank}"
        pathlib.Path(my_claim_file).touch()

        # Wait a second to make sure we don't have a collision where more than
        # one worker claimed the same directory
        time.sleep(1)

        # Make sure there is still only one claim file in this directory. If
        # more than one claim file, the file with lowest rank wins.
        claim_files = glob.glob(f"{claim_base}*")
        if len(claim_files) > 1:

            # Get ranks for all claim files.
            ranks = [int(c.split("-")[-1]) for c in claim_files]
            ranks.sort()

            # If another file exists created by a worker with a lower rank,
            # that one will do the calculation in this directory.
            if ranks[0] != rank:
                os.remove(my_claim_file)
                os.chdir("..")
                continue

        # If we got all the way here, this directory is ours to run in. Run the
        # run_generax.sh script
        subprocess.run(["bash","run_generax.sh"])

        # Job is done. (Whohoo!). Create a "completed" file.
        my_completed_file = f"{complete_base}_run-{run_id}_by-{rank}"
        pathlib.Path(my_completed_file).touch()

        # Nuke claim file and move on to next directory
        os.remove(my_claim_file)
        os.chdir("..")

    os.chdir("..")

if __name__ == "__main__":
    main()
