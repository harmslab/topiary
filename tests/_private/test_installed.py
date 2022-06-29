
import pytest

from topiary._private import installed

import warnings

def test__version_checker():

    # Good check -- should work
    cmd = ["git","--version"]
    def _version_slicer(ret):
        return ret.stdout.decode().split()[2].strip()

    v = installed._version_checker(cmd,_version_slicer)
    assert type(v) is tuple
    assert len(v) > 1
    assert int(v[0]) > 0

    # Should fail with binary not found (-2,-2,-2)
    cmd = ["not_really_a_binary"]
    def _version_slicer(ret):
        return ret.stdout.decode().split()[2].strip()

    v = installed._version_checker(cmd,_version_slicer)
    assert type(v) is tuple
    assert v == (-2,-2,-2)

    # Should fail with could not run (-1,-1,-1)
    cmd = ["git","--bad_argument_not_recognized"]
    def _version_slicer(ret):
        return ret.stdout.decode().split()[2].strip()

    v = installed._version_checker(cmd,_version_slicer)
    assert type(v) is tuple
    assert v == (-1,-1,-1)

    # Should fail with ran but could not get version (0,0,0)
    cmd = ["git","--version"]
    def _version_slicer(ret):
        # bad parsing call -- last split()[1] will throw an IndexError
        return ret.stdout.decode().split()[2].strip().split()[1]

    v = installed._version_checker(cmd,_version_slicer)
    assert type(v) is tuple
    assert v == (0,0,0)

def test_check_muscle():

    version = installed.check_muscle()

    if version == (-2,-2,-2):
        warnings.warn("muscle not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("muscle is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("muscle is installed but we cannot parse its version string!")


def test_check_generax():

    version = installed.check_generax()

    if version == (-2,-2,-2):
        warnings.warn("generax not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("generax is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("generax is installed but we cannot parse its version string!")

def test_check_raxml():

    version = installed.check_raxml()

    if version == (-2,-2,-2):
        warnings.warn("raxml-ng not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("raxml-ng is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("raxml-ng is installed but we cannot parse its version string!")

def test_check_blastp():

    version = installed.check_blastp()

    if version == (-2,-2,-2):
        warnings.warn("blastp not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("blastp is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("blastp is installed but we cannot parse its version string!")

def test_check_makeblastdb():

    version = installed.check_makeblastdb()

    if version == (-2,-2,-2):
        warnings.warn("makeblastdb not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("makeblastdb is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("makeblastdb is installed but we cannot parse its version string!")
