"""
Check for installed external software in the path.
"""

import topiary
import numpy as np
import subprocess, shutil, os, re


def _version_checker(cmd,version_slicer):
    """
    Generic version checking function.

    Parameters
    ----------
    cmd : list
        command to pass via subprocess.run
    version_slicer : function
        function that takes the subprocess return value and pulls a version
        string out of the captured stdout/stderr.

    Returns
    -------
    binary_path : str or None
        path to binary. If program is not in the path, return None
    version : tuple
        output meanings:
        + :code:`(-2,-2,-2)`, not found
        + :code:`(-1,-1,-1)`, found but does not run
        + :code:`(0,0,0)` found but could not figure out version
        + :code:`(major,minor,patch)` i.e. (3.8.1). This is done by splitting on
          the :code:`.` character, so this will always be a tuple but may have
          any length > 1. Also, the elements will be :code:`str` not :code:`int`.
    """

    # Get path to binary file
    binary_path = shutil.which(cmd[0])
    if binary_path is None:
        return None, (-2,-2,-2)

    # Run and attempt to get version
    ret = subprocess.run(cmd,capture_output=True)
    if ret.returncode != 0:
        return binary_path, (-1,-1,-1)

    try:
        version = version_slicer(ret)
        if version[0] in ["v","V"]:
            version = version[1:]
        version = tuple(version.split("."))
    except:
        return binary_path, (0,0,0)

    return binary_path, version


def check_muscle(binary=None):
    """
    Check for muscle in the PATH and get its version.

    Returns
    -------
    binary_path : str or None
        path to binary. If program is not in the path, return None
    version : tuple
        output meanings:
        + :code:`(-2,-2,-2)`, not found
        + :code:`(-1,-1,-1)`, found but does not run
        + :code:`(0,0,0)` found but could not figure out version
        + :code:`(major,minor,patch)` i.e. (3.8.1). This is done by splitting on
          the :code:`.` character, so this will always be a tuple but may have
          any length > 1. Also, the elements will be :code:`str` not :code:`int`.
    """

    def _version_slicer(ret):
        return ret.stdout.split()[1].decode()

    if binary is None:
        binary = "muscle"

    return _version_checker([binary],_version_slicer)


def check_generax(binary=None):
    """
    Check for generax in the PATH and get its version.

    Returns
    -------
    binary_path : str or None
        path to binary. If program is not in the path, return None
    version : tuple
        output meanings:
        + :code:`(-2,-2,-2)`, not found
        + :code:`(-1,-1,-1)`, found but does not run
        + :code:`(0,0,0)` found but could not figure out version
        + :code:`(major,minor,patch)` i.e. (3.8.1). This is done by splitting on
          the :code:`.` character, so this will always be a tuple but may have
          any length > 1. Also, the elements will be :code:`str` not :code:`int`.
    """

    def _version_slicer(ret):
        lines = ret.stdout.decode().split("\n")
        for line in lines:
            if re.search("generax",line,flags=re.IGNORECASE):
                return line.split()[2:][0].strip()
        return None

    if binary is None:
        binary = "generax"

    return _version_checker([binary],_version_slicer)


def check_raxml(binary=None):
    """
    Check for raxml-ng in the PATH and get its version.

    Returns
    -------
    binary_path : str or None
        path to binary. If program is not in the path, return None
    version : tuple
        output meanings:
        + :code:`(-2,-2,-2)`, not found
        + :code:`(-1,-1,-1)`, found but does not run
        + :code:`(0,0,0)` found but could not figure out version
        + :code:`(major,minor,patch)` i.e. (3.8.1). This is done by splitting on
          the :code:`.` character, so this will always be a tuple but may have
          any length > 1. Also, the elements will be :code:`str` not :code:`int`.
    """

    def _version_slicer(ret):
        return ret.stdout.decode().strip().split("\n")[0].split()[2]

    if binary is None:
        binary = "raxml-ng"

    return _version_checker([binary],_version_slicer)


def check_blastp(binary=None):
    """
    Check for blastp in the PATH and get its version.

    Returns
    -------
    binary_path : str or None
        path to binary. If program is not in the path, return None
    version : tuple
        output meanings:
        + :code:`(-2,-2,-2)`, not found
        + :code:`(-1,-1,-1)`, found but does not run
        + :code:`(0,0,0)` found but could not figure out version
        + :code:`(major,minor,patch)` i.e. (3.8.1). This is done by splitting on
          the :code:`.` character, so this will always be a tuple but may have
          any length > 1. Also, the elements will be :code:`str` not :code:`int`.
    """

    def _version_slicer(ret):
        return ret.stdout.decode().split()[1].strip()

    if binary is None:
        binary = "blastp"

    return _version_checker([binary,"-version"],_version_slicer)

def check_makeblastdb(binary=None):
    """
    Check for makeblastdb in the PATH and get its version.

    Returns
    -------
    binary_path : str or None
        path to binary. If program is not in the path, return None
    version : tuple
        output meanings:
        + :code:`(-2,-2,-2)`, not found
        + :code:`(-1,-1,-1)`, found but does not run
        + :code:`(0,0,0)` found but could not figure out version
        + :code:`(major,minor,patch)` i.e. (3.8.1). This is done by splitting on
          the :code:`.` character, so this will always be a tuple but may have
          any length > 1. Also, the elements will be :code:`str` not :code:`int`.
    """

    def _version_slicer(ret):
        return ret.stdout.decode().split()[1].strip()

    if binary is None:
        binary = "makeblastdb"

    return _version_checker([binary,"-version"],_version_slicer)

def check_git(binary=None):
    """
    Check for git in the PATH and get its version.

    Returns
    -------
    binary_path : str or None
        path to binary. If program is not in the path, return None
    version : tuple
        output meanings:
        + :code:`(-2,-2,-2)`, not found
        + :code:`(-1,-1,-1)`, found but does not run
        + :code:`(0,0,0)` found but could not figure out version
        + :code:`(major,minor,patch)` i.e. (3.8.1). This is done by splitting on
          the :code:`.` character, so this will always be a tuple but may have
          any length > 1. Also, the elements will be :code:`str` not :code:`int`.
    """

    def _version_slicer(ret):
        return ret.stdout.decode().split()[2].strip()

    if binary is None:
        binary = "git"

    return _version_checker([binary,"--version"],_version_slicer)

def check_mpirun(binary=None):
    """
    Check for mpirun in the PATH and get its version.

    Returns
    -------
    binary_path : str or None
        path to binary. If program is not in the path, return None
    version : tuple
        output meanings:
        + :code:`(-2,-2,-2)`, not found
        + :code:`(-1,-1,-1)`, found but does not run
        + :code:`(0,0,0)` found but could not figure out version
        + :code:`(major,minor,patch)` i.e. (3.8.1). This is done by splitting on
          the :code:`.` character, so this will always be a tuple but may have
          any length > 1. Also, the elements will be :code:`str` not :code:`int`.
    """

    def _version_slicer(ret):
        return ret.stdout.decode().split("\n")[0].split()[-1]

    if binary is None:
        binary = "mpirun"

    return _version_checker([binary,"--version"],_version_slicer)

def _compare_versions(installed,required):
    """
    Compare an installed version tuple to a required version tuple.

    Parameters
    ----------
    installed : tuple
        tuple returned by check_xxx. elements are str.
    required : tuple
        required version as a tuple. elements are int. examples: (1,1) would
        require version 1.1. (1,) would require at least version 1.

    Returns
    -------
    status : bool or None
        True (version is high enough); False (version is not high enough);
        None (not clear if version is high enough or not)
    """

    # Convert up to len(required) elements of installed version
    # into integers for comparison
    tmp = []
    length = 0
    for i in range(len(required)):
        try:
            tmp.append(int(installed[i]))
            length += 1
        except (TypeError,IndexError,ValueError):
            break

    # Figure out if we compared all elements specified in the
    # required tuple
    full_length = False
    if length == len(required):
        full_length = True

    # Compare shared elements between installed and required
    installed = np.array(tmp)[:length]
    required = np.array(required)[:length]
    difference = installed - required

    # Go through elements
    for d in difference:

        # Installed better than required -- return succcess
        if d > 0:
            return True

        # Installed same as required -- move to next version level
        if d == 0:
            continue

        # Installed worse than required -- return failure
        if d < 0:
            return False

    # If we could parse full length of installed version and we got here,
    # version matches exactly.
    if full_length:
        return True

    # If we got here, we did not have the full length of required version
    # but those that were present matched. Ambiguous.
    return None

def validate_stack(to_check):

    binary_tests = {"makeblastdb":check_makeblastdb,
                    "blastp":check_blastp,
                    "raxml-ng":check_raxml,
                    "generax":check_generax,
                    "muscle":check_muscle,
                    "git":check_git,
                    "mpirun":check_mpirun}

    out = []
    bad_prog = []
    for check in to_check:

        program = check["program"]
        min_version = check["min_version"]
        must_pass = check["must_pass"]

        try:
            binary = check["binary"]
        except KeyError:
            binary = None

        out.append(70*"-")
        out.append(f"Checking {binary}")
        out.append(70*"-")
        out.append("")

        fcn = binary_tests[program]

        binary, version = fcn(binary)

        if version == (-2,-2,-2):
            installed =   "N"
            binary_path = "-"
            binary_runs = "-"
            version_str = "-"
            passes =      "N"

        elif version == (-1,-1,-1):
            installed =   "Y"
            binary_path = binary
            binary_runs = "N"
            version_str = "-"
            passes =      "N"

        elif version == (0,0,0):
            installed =   "Y"
            binary_path = binary
            binary_runs = "Y"
            version_str = "-"
            passes =      "?"

        else:
            installed =   "Y"
            binary_path = binary
            binary_runs = "Y"
            version_str = ".".join(version)

            status = _compare_versions(version,min_version)
            if status is True:
                passes = "Y"
            elif status is False:
                passes = "N"
            else:
                passes = "?"

        min_version_str = ".".join([str(v) for v in min_version])

        out.append(f"    installed:       {installed}")
        out.append(f"    binary_path:     {binary_path}")
        out.append(f"    binary runs:     {binary_runs}")
        out.append(f"    version:         {version_str}")
        out.append(f"    minimum version: {min_version_str}")
        out.append(f"    passes:          {passes}")

        if passes == "N" or (passes == "?" and must_pass):
            bad_prog.append((program,min_version_str))

        out.append("")

    print("\n".join(out),flush=True)

    if len(bad_prog) > 0:
        err = "\nNot all programs available. Please make sure that the following\n"
        err += "programs are in the $PATH.\n"
        for b, v in bad_prog:
            err += f" + {b}>={v}\n"
        err += "\n"
        err += "The current $PATH visible to python is:\n"
        err += f"    {os.environ['PATH']}"
        err += "\n"
        raise RuntimeError(err)

def test_mpi_configuration(num_threads):
    """
    Make sure mpi configuration allows the requested number of threads.
    """

    # if threads were not passed in directly, infer from the environment
    if num_threads == -1:
        num_threads = topiary._private.threads.get_num_threads(num_threads)

    # Run ls on num_threads.
    cmd = ["mpirun","-np",f"{num_threads}","ls"]
    ret = subprocess.run(cmd,capture_output=True)

    # If mpirun failed,
    if ret.returncode != 0:

        err = "mpirun is not working. See error below. This could because you\n"
        err += "set --num_threads to be more than the number of nodes you have\n"
        err += "allocated on your cluster. If you did not set --num_threads\n"
        err += "specifically, try setting it rather than having topiary try to\n"
        err += "figure out the number of processors. Another issue could be subtle\n"
        err += "problems with how processors are being requested via your job\n"
        err += "management software. Maybe play with flags like --ntasks-per-node\n"
        err += "or talk to your cluster administrator. mpirun stdout and stderr\n"
        err += "follows:\n\n"
        err += "stdout:\n\n"
        err += f"{ret.stdout.decode()}"
        err += "\nstderr:\n\n"
        err += f"{ret.stderr.decode()}"
        err += "\n"

        raise ValueError(err)
