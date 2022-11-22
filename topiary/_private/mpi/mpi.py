"""
Functions for interfacing with and validating mpi configuration.
"""

import subprocess
import sys
import os

from topiary._private.environment import load_env_variable
from topiary._private.check import check_int

def get_hosts(num_slots):
    """
    Get the hosts allocated by the job manager by running a script with mpirun
    that returns the host names.

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
    if ret.returncode == 0:
        stdout = ret.stdout.decode()
        hosts = [s.strip() for s in stdout.split("\n") if s.strip() != ""]
        hosts.sort()
    else:
        err = "Could not determine hosts. _get_hosts.py script returned:\n\n"
        err += f"stdout:\n\n{ret.stdout.decode()}\n\n"
        err += f"stderr:\n\n{ret.stderr.decode()}\n\n"
        raise RuntimeError(err)

    return hosts

def get_num_slots():
    """
    Get the number of mpi slots available by running get_hosts until it throws
    an error. 

    Returns
    -------
    num_slots : int
        number of slots available
    """

    # Get environment variable if defined
    max_num_slots = load_env_variable("TOPIARY_MAX_SLOTS",
                                      check_function=check_int,
                                      check_function_kwargs={"minimum_allowed":1})


    # Increase number of slots until mpirun fails
    num_slots = 1
    while True:

        try:
            _ = get_hosts(num_slots)
            num_slots += 1
        except RuntimeError as error:
            num_slots = num_slots - 1

            # If we have no slots, something went terribly wrong.
            if num_slots < 1:
                err = "\nCould not determine the number of MPI slots. The test script\n"
                err += "raised the following error."
                raise RuntimeError(err) from error

            break

        # Hard cap based on environment variable
        if max_num_slots is not None:
            if num_slots >= max_num_slots:
                num_slots = max_num_slots
                break


    return num_slots


def check_mpi_configuration(num_threads):
    """
    Make sure mpi configuration allows the requested number of threads.

    Parameters
    ----------
    num_threads : int
        number of threads (e.g. slots) to test. if -1, try to infer the number
        of slots using get_num_slots
    """

    # if threads were not passed in directly, infer from the environment
    if num_threads == -1:
        num_threads = get_num_slots()

    try:
        get_hosts(num_threads)
    except RuntimeError as error:
        err = "\n\nmpirun is not working. This could because you\n"
        err += "set num_threads to be more than the number of nodes you have\n"
        err += "allocated. If you did not set num_threads specifically, try\n"
        err += "setting it rather than having topiary try to figure out the\n"
        err += "number of processors. Another issue could be subtle problems\n"
        err += "with how processors are being requested via your job management\n"
        err += "software (i.e. SLURM, TORQUE, etc.). Maybe play with flags like\n"
        err += "--ntasks-per-node or talk to your cluster administrator.\n"

        raise RuntimeError(err) from error
