"""
Functions for interfacing with and validating mpi configuration.
"""

import subprocess
import sys
import os

def get_hosts(num_slots):
    """
    Get the hosts allocated by the job manager by running a test script with
    mpirun that then returns the host names.

    Returns
    -------
    host : list
        list of host names (for example: ["n001","n001","n002"])
    """

    location = os.path.dirname(os.path.realpath(__file__))
    script = os.path.join(location,"_get_hosts.py")
    python = sys.executable


    ret = subprocess.run(["mpirun","-np",f"{num_slots}",python,script],
                          capture_output=True)
    stdout = ret.stdout.decode()
    hosts = [s.strip() for s in stdout.split("\n") if s.strip() != ""]
    hosts.sort()

    return hosts

def get_num_slots(test_binary):
    """
    Get the number of mpi slots available.

    Parameters
    ----------
    test_binary : str
        path to an mpi enabled binary file

    Returns
    -------
    num_slots : int
        number of slots available
    """

    # Increase number of slots until mpirun fails
    num_slots = 1
    while True:

        cmd = ["mpirun","-np",f"{num_slots}",test_binary]
        ret = subprocess.run(cmd,capture_output=True)
        if ret.returncode == 0:
            num_slots += 1
            continue

        break

    num_slots = num_slots - 1

    # If we have no slots, something went terribly wrong.
    if num_slots < 1:
        err = "\nCould not determine the number of MPI slots.\n"
        err += f"Tried to run: {' '.join(cmd)}\n"
        err += "This returned:\n\n"
        err += "stdout:\n\n"
        err += f"{ret.stdout.decode()}"
        err += "\nstderr:\n\n"
        err += f"{ret.stderr.decode()}"
        err += "\n"
        raise RuntimeError(err)

    return num_slots


def check_mpi_configuration(num_threads,test_binary):
    """
    Make sure mpi configuration allows the requested number of threads.

    Parameters
    ----------
    num_threads : int
        number of threads (e.g. slots) to test. if -1, try to infer the number
        of slots using get_num_slots
    test_binary : str
        path to an mpi enabled binary file
    """

    # if threads were not passed in directly, infer from the environment
    if num_threads == -1:
        num_threads = get_num_slots(test_binary)

    # Run ls on num_threads.
    cmd = ["mpirun","-np",f"{num_threads}",test_binary]
    ret = subprocess.run(cmd,capture_output=True)

    # If mpirun failed,
    if ret.returncode != 0:

        err = "\n\nmpirun is not working. See error below. This could because you\n"
        err += "set num_threads to be more than the number of nodes you have\n"
        err += "allocated. If you did not set num_threads specifically, try\n"
        err += "setting it rather than having topiary try to figure out the\n"
        err += "number of processors. Another issue could be subtle problems\n"
        err += "with how processors are being requested via your job management\n"
        err += "software (i.e. SLURM, TORQUE, etc.). Maybe play with flags like\n"
        err += "--ntasks-per-node or talk to your cluster administrator. mpirun\n"
        err += "stdout and stderr follows:\n\n"
        err += "stdout:\n\n"
        err += f"{ret.stdout.decode()}"
        err += "\nstderr:\n\n"
        err += f"{ret.stderr.decode()}"
        err += "\n"

        raise ValueError(err)
