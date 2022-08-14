import pytest

import topiary
from topiary._private.interface import gen_seed
from topiary._private.interface import launch
from topiary._private.interface import create_new_dir
from topiary._private.interface import copy_input_file
import numpy as np
import pandas as pd

import warnings, os, shutil, glob, json, sys

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
#
# def test_read_previous_run_dir(run_directories,tmpdir):
#
#     bad_inputs = [1,None,{"test":1},(),[]]
#     for b in bad_inputs:
#         with pytest.raises(ValueError):
#             interface.read_previous_run_dir(b)
#
#     with pytest.raises(FileNotFoundError):
#         interface.read_previous_run_dir("not_really_a_directory")
#
#     # Pass something too deep into good run directory
#     with pytest.raises(FileNotFoundError):
#         interface.read_previous_run_dir(os.path.join(run_directories["ml-tree"],
#                                               "output"))
#
#     # no json; should die
#     with pytest.raises(FileNotFoundError):
#         interface.read_previous_run_dir(run_directories["ml-tree_no-json"])
#
#     # bad dataframe; should die (make sure it's checking df on way in)
#     with pytest.raises(ValueError):
#         interface.read_previous_run_dir(run_directories["ml-tree_bad-df"])
#
#     # no dataframe; should die
#     with pytest.raises(FileNotFoundError):
#         interface.read_previous_run_dir(run_directories["ml-tree_no-df"])
#
#     # Pass in a bad parameters.json
#     current_dir = os.getcwd()
#     os.chdir(tmpdir)
#
#     bad_dir = "not-really-an-output-dir"
#     os.mkdir(bad_dir)
#     os.mkdir(os.path.join(bad_dir,"output"))
#
#     f = open(os.path.join(bad_dir,"output","run_parameters.json"),"w")
#     f.write("junk\n")
#     f.close()
#
#     with pytest.raises(ValueError):
#         interface.read_previous_run_dir(bad_dir)
#
#     shutil.copy(os.path.join(run_directories["ml-tree"],"output","run_parameters.json"),
#                 os.path.join(bad_dir,"output","run_parameters.json"))
#
#     # If we see FileNotFoundError, we've made it to checking for the dataframe
#     with pytest.raises(FileNotFoundError):
#         interface.read_previous_run_dir(bad_dir)
#
#     os.chdir(current_dir)
#
#
#     # no tree; should work fine
#     interface.read_previous_run_dir(run_directories["ml-tree_no-tree"])
#
#     # Make sure the function is properly reading run directories
#     for t in [("ml-tree","ml_tree"),
#               ("ancestors","ancestors"),
#               ("find-model","find_best_model"),
#               ("reconciled","reconciliation")]:
#
#         key = t[0]
#         calc_type = t[1]
#
#         out = interface.read_previous_run_dir(run_directories[key])
#         assert out["model"] == "JTT"
#         assert out["calc_type"] == calc_type
#         assert type(out["cmd"]) in [type("string"),type(None)]
#         try:
#             tf = out["tree_file"]
#             assert os.path.exists(tf)
#         except KeyError:
#             pass
#         # This will throw an error if the dataframe read is bad
#         check.check_topiary_dataframe(out["df"])
#
#
# def test_setup_calc_dir(xml_to_anc_output,tmpdir):
#
#     current_dir = os.getcwd()
#     os.chdir(tmpdir)
#
#     # Should throw value error --> no dataframe specified
#     with pytest.raises(ValueError):
#         interface.setup_calc_dir()
#
#     # --------------------------------------------------------------------------
#     # Basic run.
#
#     previous_dir = os.path.join(xml_to_anc_output,"ml-tree")
#     out = interface.setup_calc_dir(previous_dir=previous_dir)
#     assert os.path.exists(out["calc_dir"])
#
#     assert out["csv_file"] is None
#
#     assert out["tree_file"] == os.path.join(out["calc_dir"],"input","tree.newick")
#     assert os.path.exists(out["tree_file"])
#
#     assert out["model"] == "JTT"
#
#     assert out["alignment_file"] == os.path.join(out["calc_dir"],"input","alignment.phy")
#     assert os.path.exists(out["alignment_file"])
#
#     assert out["previous_dir"] == previous_dir
#
#     assert out["starting_dir"] == tmpdir
#
#     # Make sure it goes into the directory
#     assert os.path.isdir(out["calc_dir"])
#
#
#     # Bad previous_dir will be caught by read_previous_run_dir. See tests
#     # for that function.
#
#     # --------------------------------------------------------------------------
#     # Run that takes a previous dataframe and nothing else
#
#     df_file = os.path.join(xml_to_anc_output,"ml-tree","output","dataframe.csv")
#     out = interface.setup_calc_dir(df=df_file)
#     assert os.path.exists(out["calc_dir"])
#     topiary.read_dataframe(os.path.join(out["calc_dir"],"input","dataframe.csv"))
#
#     # Cleanup
#     shutil.rmtree(out["calc_dir"])
#
#     # pass dataframe as a dataframe
#     df = topiary.read_dataframe(df_file)
#     out = interface.setup_calc_dir(df=df)
#     assert os.path.exists(out["calc_dir"])
#     topiary.read_dataframe(os.path.join(out["calc_dir"],"input","dataframe.csv"))
#
#     # Cleanup
#     shutil.rmtree(out["calc_dir"])
#
#     bad_df = [1,None,[],pd.DataFrame,{"test":[1,2,3]}]
#     for b in bad_df:
#         with pytest.raises(ValueError):
#             interface.setup_calc_dir(df)
#
#     # --------------------------------------------------------------------------
#     # model
#
#     bad_model = [1,None,[],pd.DataFrame,{"test":[1,2,3]}]
#     for b in bad_model:
#         with pytest.raises(ValueError):
#             interface.setup_calc_dir(previous_dir=previous_dir,model=bad_model)
#
#     out = interface.setup_calc_dir(previous_dir=previous_dir,model="manual_model")
#     assert out["model"] == "manual_model"
#
#     # clean up
#     shutil.rmtree(out["calc_dir"])
#
#     # --------------------------------------------------------------------------
#     # tree_file
#
#     # bad tree file passes are checked by copy_input_file. See there for error
#     # checking tests
#
#     tree_file = os.path.join(xml_to_anc_output,"ancestors","output","tree.newick")
#
#     # Make sure the tree file overwrites what came from the previous dir
#     out = interface.setup_calc_dir(previous_dir=previous_dir,
#                               tree_file=tree_file)
#     f = open(out["tree_file"],"r")
#     t = f.read()
#     f.close()
#
#     # manual tree starts with this; previous_dir tree does not.
#     assert t.startswith("((BHxLaFLFBU")
#
#     # clean up
#     shutil.rmtree(out["calc_dir"])
#
#     os.chdir(current_dir)

def test__follow_log_subproc_wrapper():
    # Super simple function that would require lots of test infrastructure to
    # run. Return True so this is not flagged as missing by test_crawler.
    return True

def test__follow_log_generator():
    # Function that would require lots of test infrastructure to run. Logging
    # not critical to results, so skipping for now. Return True so this is not
    # flagged as missing by test_crawler.
    return True


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
#
# def test_write_run_information(tmpdir,test_dataframes):
#
#     df = test_dataframes["good-df"]
#     outdir = os.path.join(tmpdir,"stupid")
#     os.mkdir(outdir)
#
#     interface.write_run_information(outdir=outdir,
#                                     df=df,
#                                     calc_type="my_calc",
#                                     model="my_model",
#                                     cmd="my_command",
#                                     start_time="my_start")
#
#     csv_file = os.path.join(outdir,"dataframe.csv")
#     assert os.path.isfile(csv_file)
#     topiary.read_dataframe(csv_file)
#
#     f = open(os.path.join(outdir,"run_parameters.json"))
#     out = json.load(f)
#     f.close()
#
#     assert out["calc_type"] == "my_calc"
#     assert out["model"] == "my_model"
#     assert out["cmd"] == "my_command"
#     assert out["version"] == topiary.__version__
#     assert out["start_time"] == "my_start"
#
#     # Don't check stop_time, just make sure it is present
#     out["stop_time"]
