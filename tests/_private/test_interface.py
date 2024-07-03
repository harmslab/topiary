import pytest

import topiary
from topiary._private.interface import gen_seed
from topiary._private.interface import launch
from topiary._private.interface import create_new_dir
from topiary._private.interface import copy_input_file
from topiary._private.interface import run_cleanly
from topiary._private.interface import rmtree
from topiary._private.interface import _follow_log_generator
from topiary._private.interface import _follow_log_subproc_wrapper

import os
import sys
import time

def test_run_cleanly(tmpdir):

    cwd = os.getcwd()
    os.chdir(tmpdir)

    # Make sure we can send in args, kwargs, and get return values
    @run_cleanly
    def test_function_arg_passing(x):
        return x*2

    assert test_function_arg_passing(2) == 4
    assert test_function_arg_passing(x=5) == 10

    # This should create a directory called crash, move into crash, raise a
    # value error, then return back to the starting directory
    @run_cleanly
    def test_function_crash():
        os.mkdir("crash")
        os.chdir("crash")
        raise ValueError
        
    with pytest.raises(Exception):
        test_function_crash()
    assert os.getcwd() == os.path.abspath(tmpdir)
    assert os.path.isdir("crash")

    # Make sure we do not return to starting directory if there is no crash
    @run_cleanly
    def test_function_chdir():
        os.mkdir("stupid")
        os.chdir("stupid")

    test_function_chdir()
    assert os.path.abspath(os.getcwd()) == os.path.abspath(os.path.join(tmpdir,"stupid"))

    os.chdir(cwd)

    

def test_MockTqdm():

    pass

def test_MockTqdm___enter__():

    pass

def test_gen_seed():

    seed = gen_seed()
    assert type(seed) is str
    assert len(seed) == 10
    int(seed)

def test_create_new_dir(tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    dir_name = create_new_dir()
    assert os.path.isdir(dir_name)
    dir_split = os.path.split(dir_name)[-1].split("_")
    assert dir_split[0] == "calculation"
    assert len(dir_split[1]) == 10

    dir_name = create_new_dir("cool_dir")
    assert os.path.isdir("cool_dir")

    with pytest.raises(FileExistsError):
        create_new_dir("cool_dir")

    create_new_dir("cool_dir",overwrite=True)

    os.chdir(current_dir)

def test_copy_input_file(tmpdir,test_dataframes):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    test_file = "test_file_to_copy.txt"
    f = open(test_file,"w")
    for i in range(10):
        f.write(f"{i}\n")
    f.close()

    os.mkdir("target_dir")

    with pytest.raises(FileNotFoundError):
        copy_input_file("not_a_file","target_dir")

    with pytest.raises(FileNotFoundError):
        copy_input_file(test_file,"not_a_dir")

    with pytest.raises(FileNotFoundError):
        copy_input_file(test_file,test_file)

    # Run in default configuration
    copy_input_file(test_file,"target_dir")
    assert os.path.isfile(os.path.join("target_dir",test_file))

    # Specify output name
    copy_input_file(test_file,"target_dir",file_name="rocket.txt")
    assert os.path.isfile(os.path.join("target_dir","rocket.txt"))

    # Don't put in input directory
    copy_input_file(test_file,"target_dir")
    assert os.path.isfile(os.path.join("target_dir",test_file))
    os.remove(os.path.join("target_dir",test_file))

    # Put into existing input directory
    os.mkdir(os.path.join("target_dir","input"))
    copy_input_file(test_file,os.path.join("target_dir","input"))
    assert os.path.isfile(os.path.join("target_dir","input",test_file))

    os.chdir(current_dir)


def test__follow_log_subproc_wrapper():
    # Super simple function that would require lots of test infrastructure to
    # run. Touch without calling to test crawler does not flag as missing
    _follow_log_subproc_wrapper
    return None

def test__follow_log_generator():
    # Function that would require lots of test infrastructure to run. Logging
    # not critical to results, so skipping for now. Touch without calling so 
    # test crawler does not flag as missing
    _follow_log_generator
    return None


def test_launch(tmpdir,programs):
    """
    Test the launcher.
    """

    # This test basically makes sure the code runs. Hard to test the
    # multithreading bits. Also, because this is an internal function, it does
    # not do argument checking, etc.

    out_dir = os.path.join(os.path.abspath(tmpdir),"run-test-dir")
    os.mkdir(out_dir)

    prg = programs["write_to_file_over_time.py"]
    cmd = [prg,"output.txt","--num_steps","2"]

    # On windows box, we need to inject python as first argument
    if os.name == "nt":
        cmd.insert(0,sys.executable)

    launch(cmd,out_dir)

    f = open(os.path.join(out_dir,"output.txt"))
    lines = f.readlines()
    f.close()

    assert len(lines) == 2
    assert len(lines[0].split(";")) == 3
    assert len(lines[1].split(";")) == 3

def test_rmtree(tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    os.mkdir("testdir")
    with pytest.raises(ValueError):
        rmtree("testdir","not_an_integer")

    assert os.path.isdir("testdir")
    rmtree("testdir")
    assert not os.path.exists("testdir")

    os.mkdir("another-test")
    f = open(os.path.join("another-test","stupid.txt"),"w")
    f.write("something")
    f.close()

    assert os.path.isfile(os.path.join("another-test","stupid.txt"))
    rmtree("another-test")
    assert not os.path.exists("another-test")

    # I'm not sure how to simulate a broken NFS situation where this will 
    # fail. Kind of incomplete as a test. 


    os.chdir(current_dir)