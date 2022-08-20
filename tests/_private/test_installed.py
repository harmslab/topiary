
import pytest

from topiary._private.installed import _version_checker
from topiary._private.installed import check_git
from topiary._private.installed import check_muscle
from topiary._private.installed import check_generax
from topiary._private.installed import check_raxml
from topiary._private.installed import check_blastp
from topiary._private.installed import check_makeblastdb
from topiary._private.installed import check_mpirun
from topiary._private.installed import _compare_versions
from topiary._private.installed import validate_stack

from topiary.generax import GENERAX_BINARY

import warnings
import os

def test__version_checker():

    # Good check -- should work
    cmd = ["git","--version"]
    def _version_slicer(ret):
        return ret.stdout.decode().split()[2].strip()

    b, v = _version_checker(cmd,_version_slicer)
    assert type(b) is str
    assert type(v) is tuple
    assert len(v) > 1
    assert int(v[0]) > 0

    # Should fail with binary not found (-2,-2,-2)
    cmd = ["not_really_a_binary"]
    def _version_slicer(ret):
        return ret.stdout.decode().split()[2].strip()

    b, v = _version_checker(cmd,_version_slicer)
    assert b is None
    assert type(v) is tuple
    assert v == (-2,-2,-2)

    # Should fail with could not run (-1,-1,-1)
    cmd = ["git","--bad_argument_not_recognized"]
    def _version_slicer(ret):
        return ret.stdout.decode().split()[2].strip()

    b, v = _version_checker(cmd,_version_slicer)
    assert type(b) is str
    assert type(v) is tuple
    assert v == (-1,-1,-1)

    # Should fail with ran but could not get version (0,0,0)
    cmd = ["git","--version"]
    def _version_slicer(ret):
        # bad parsing call -- last split()[1] will throw an IndexError
        return ret.stdout.decode().split()[2].strip().split()[1]

    b, v = _version_checker(cmd,_version_slicer)
    assert type(b) is str
    assert type(v) is tuple
    assert v == (0,0,0)

def test_check_muscle():

    binary, version = check_muscle()

    if version == (-2,-2,-2):
        warnings.warn("muscle not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("muscle is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("muscle is installed but we cannot parse its version string!")

@pytest.mark.run_generax
def test_check_generax():

    binary, version = check_generax()

    if version == (-2,-2,-2):
        warnings.warn("generax not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("generax is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("generax is installed but we cannot parse its version string!")

@pytest.mark.run_raxml
def test_check_raxml():

    binary, version = check_raxml()

    if version == (-2,-2,-2):
        warnings.warn("raxml-ng not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("raxml-ng is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("raxml-ng is installed but we cannot parse its version string!")

def test_check_blastp():

    binary, version = check_blastp()

    if version == (-2,-2,-2):
        warnings.warn("blastp not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("blastp is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("blastp is installed but we cannot parse its version string!")

def test_check_makeblastdb():

    binary, version = check_makeblastdb()

    if version == (-2,-2,-2):
        warnings.warn("makeblastdb not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("makeblastdb is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("makeblastdb is installed but we cannot parse its version string!")

def test_check_git():

    binary, version = check_git()

    if version == (-2,-2,-2):
        warnings.warn("git not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("git is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("git is installed but we cannot parse its version string!")

# mark for run_generax because that's why we worry about MPI in the first place.
# this will avoid a warning on windows tests. 
@pytest.mark.run_generax
def test_check_mpirun():

    binary, version = check_mpirun()

    if version == (-2,-2,-2):
        warnings.warn("mpirun not installed -- skipping test")

    if version == (-1,-1,-1):
        raise RuntimeError("mpirun is installed but not working!")

    if version == (0,0,0):
        raise RuntimeError("mpirun is installed but we cannot parse its version string!")


def test__compare_versions():

    # good version string; only specified element (1,) matches
    out = _compare_versions(("1","0"),(1,))
    assert out is True

    # bad version string; only specified element (1,) matches
    out = _compare_versions(("1","1b"),(1,))
    assert out is True

    # can't compare last element -- ambiguous
    out = _compare_versions(("1","1b"),(1,1))
    assert out is None

    # first comparable element bad -- should not pass
    out = _compare_versions(("1","1b"),(2,1))
    assert out is False

    # good version string; second element too low
    out = _compare_versions(("1","0"),(1,2))
    assert out is False

    # bad version string; can't check second position
    out = _compare_versions(("1","1b"),(1,2))
    assert out is None

    # First position too low
    out = _compare_versions(("0","1b"),(1,2))
    assert out is False

    # Shoudl pass -- last element high enough
    out = _compare_versions(("1","1","1"),(1,1,0))
    assert out is True

    # Should pass, matches version exactly
    out = _compare_versions(("1","1","1"),(1,1,1))
    assert out is True


def test_validate_stack():

    # not an amazing test, but at least checks core logic of whether or not
    # version is high enough.

    validate_stack([{"program":"git",
                               "min_version":(0,0,1),
                               "must_pass":True}])

    with pytest.raises(RuntimeError):
        validate_stack([{"program":"git",
                                   "min_version":(10000000,0,1),
                                   "must_pass":True}])
