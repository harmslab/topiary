"""
Check for installed external software in the path.
"""

import subprocess

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
    version : tuple
        output meanings:
        + :code:`(-2,-2,-2)`, not found
        + :code:`(-1,-1,-1)`, found but does not run
        + :code:`(0,0,0)` found but could not figure out version
        + :code:`(major,minor,patch)` i.e. (3.8.1). This is done by splitting on
          the :code:`.` character, so this will always be a tuple but may have
          any length > 1. Also, the elements will be :code:`str` not :code:`int`.
    """

    try:
        ret = subprocess.run(cmd,capture_output=True)
    except FileNotFoundError:
        return (-2,-2,-2)

    if ret.returncode != 0:
        return (-1,-1,-1)

    try:
        version = version_slicer(ret)
        if version[0] in ["v","V"]:
            version = version[1:]
        version = tuple(version.split("."))
    except:
        return (0,0,0)

    return version


def check_muscle():
    """
    Check for muscle in the PATH and get its version.

    Returns
    -------
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

    return _version_checker(["muscle"],_version_slicer)


def check_generax():
    """
    Check for generax in the PATH and get its version.

    Returns
    -------
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
        return ret.stdout.decode().split("\n")[0].split()[2:][0]

    return _version_checker(["generax"],_version_slicer)


def check_raxml():
    """
    Check for raxml-ng in the PATH and get its version.

    Returns
    -------
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

    return _version_checker(["raxml-ng"],_version_slicer)


def check_blastp():
    """
    Check for blastp in the PATH and get its version.

    Returns
    -------
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

    return _version_checker(["blastp","-version"],_version_slicer)

def check_makeblastdb():
    """
    Check for makeblastdb in the PATH and get its version.

    Returns
    -------
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

    return _version_checker(["makeblastdb","-version"],_version_slicer)
