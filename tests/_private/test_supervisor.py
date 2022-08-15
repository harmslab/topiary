
import pytest

import topiary
from topiary._private.supervisor import Supervisor

import pandas as pd

import os
import time
import json
import shutil
import glob
import pathlib
import ete3

def test_Supervisor(test_dataframes,tmpdir,generax_data):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    df_csv = test_dataframes["good-df_real-alignment"]
    gene_tree = os.path.join(generax_data["prev-ml-run"],"01_ml-tree",
                             "output","tree.newick")

    # Not really a species tree, but fine for copying newick in and out
    species_tree = os.path.join(generax_data["prev-ml-run"],"01_ml-tree",
                                "output","tree.newick")

    # Not really a reconciled tree, but fine for copying newick in and out
    reconciled_tree = os.path.join(generax_data["prev-ml-run"],"01_ml-tree",
                                "output","tree.newick")


    sv = Supervisor()

    assert sv.run_parameters["version"] == topiary.__version__
    assert sv.status == "empty"
    assert sv.run_parameters["calc_status"] == "empty"

    assert sv.input_dir is None
    assert sv.working_dir is None
    assert sv.output_dir is None
    assert sv.starting_dir == tmpdir
    assert sv.gene_tree is None
    assert sv.species_tree is None
    assert sv.reconciled_tree is None
    assert sv.alignment is None
    assert isinstance(sv.seed,int)
    assert sv.run_parameters["seed"] == sv.seed
    assert sv.df is None
    assert sv.previous_entries is None
    assert len(sv.run_parameters) == 4

    # Test different ways to generate seeds
    sv = Supervisor()
    assert isinstance(sv.seed,int)

    sv = Supervisor(seed=12345)
    assert sv.seed == 12345

    with pytest.raises(ValueError):
        Supervisor(seed="stpuid")

    # Create non-trivial past directory so we can pass in as calc_dir
    sv = Supervisor()
    sv.create_calc_dir("test0",
                       calc_type="something",
                       df=df_csv,
                       gene_tree=gene_tree)
    shutil.copy(os.path.join(sv.input_dir,"dataframe.csv"),
                os.path.join(sv.output_dir,"dataframe.csv"))
    shutil.copy(os.path.join(sv.input_dir,"gene-tree.newick"),
                os.path.join(sv.output_dir,"gene-tree.newick"))
    sv.event("test")
    sv.update("some_entry","another")
    sv.update("model","another")
    sv.finalize()

    # Get model and seed from this pass
    old_model = sv.model
    old_seed = sv.seed

    # Load that previous directory
    sv = Supervisor("test0")
    assert sv.model == old_model
    assert sv.seed == old_seed
    assert os.path.isfile(os.path.join(sv.output_dir,"dataframe.csv"))
    assert os.path.isfile(os.path.join(sv.output_dir,"gene-tree.newick"))
    assert os.path.isfile(os.path.join(sv.input_dir,"dataframe.csv"))
    assert os.path.isfile(os.path.join(sv.input_dir,"gene-tree.newick"))
    assert sv.run_parameters["calc_type"] == "something"
    assert sv.run_parameters["some_entry"] == "another"
    assert sv.previous_entries is None
    assert len(sv.run_parameters["events"]) == 1

    os.chdir(current_dir)

def test_Supervisor__load_existing(tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)
    sv = Supervisor()

    os.mkdir("junk")
    sv._calc_dir = "junk"

    with pytest.raises(FileNotFoundError):
        sv._load_existing()

    f = open(os.path.join("junk","run_parameters.json"),"w")
    f.write("yo")
    f.close()
    with pytest.raises(ValueError):
        sv._load_existing()

    f = open(os.path.join("junk","run_parameters.json"),"w")
    json.dump({"test":1},f)
    f.close()

    sv._load_existing()

    os.chdir(current_dir)

def test_Supervisor__increment(tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    sv.create_calc_dir("test0","test_no_prev")
    assert sv.previous_entries is None
    this_seed = sv.seed
    sv.event('somewhere')
    sv.update("model","mock")
    sv.update("something","else")
    sv.finalize()

    sv.create_calc_dir("test1","test_no_prev")

    # These should have been preserved
    assert sv.run_parameters["model"] == "mock"
    assert sv.run_parameters["version"] == topiary.__version__
    assert sv.run_parameters["seed"] == this_seed

    # Should *not* have been preserved
    assert "something" not in sv.run_parameters
    assert len(sv.previous_entries) == 1
    assert sv.previous_entries[0]["something"] == "else"

    os.chdir(current_dir)

def test_Supervisor_create_calc_dir(test_dataframes,tmpdir,generax_data):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    df_csv = test_dataframes["good-df_real-alignment"]
    gene_tree = os.path.join(generax_data["prev-ml-run"],"01_ml-tree",
                             "output","tree.newick")

    # Not really a species tree, but fine for copying newick in and out
    species_tree = os.path.join(generax_data["prev-ml-run"],"01_ml-tree",
                                "output","tree.newick")

    # Not really a reconciled tree, but fine for copying newick in and out
    reconciled_tree = os.path.join(generax_data["prev-ml-run"],"01_ml-tree",
                                "output","tree.newick")

    sv = Supervisor()
    assert sv.status == "empty"

    sv.create_calc_dir("test0x","test_no_prev")
    assert os.path.isdir(sv.input_dir)
    assert os.path.isdir(sv.working_dir)
    assert os.path.isdir(sv.output_dir)
    assert os.path.isfile(os.path.join("test0x","run_parameters.json"))
    assert sv.gene_tree is None
    assert sv.species_tree is None
    assert sv.reconciled_tree is None
    assert sv.alignment is None
    assert sv.starting_dir == tmpdir
    assert sv.run_parameters["calc_type"] == "test_no_prev"
    assert sv.status == "running"
    assert sv.run_parameters["creation_time"] <= time.time()

    json_file = os.path.join(sv.calc_dir,"run_parameters.json")
    f = open(json_file)
    p = json.load(f)
    f.close()

    assert p["calc_status"] == "running"
    assert p["seed"] == sv.seed
    assert p["calc_type"] == "test_no_prev"

    # should not work because dir already exists
    sv = Supervisor()
    with pytest.raises(FileExistsError):
        sv.create_calc_dir("test0x","test_no_prev_2",df=df_csv)

    # overwrite
    sv = Supervisor()
    sv.create_calc_dir("test0x","test_no_prev_2",overwrite=True,df=df_csv)
    assert sv.run_parameters["calc_type"] == "test_no_prev_2"

    # Declare it had an error, then try to make new. should throw error we can
    # overcome via force
    sv.finalize(successful=False)
    with pytest.raises(ValueError):
        sv.create_calc_dir("test1x","test_no_prev_3")
    sv.create_calc_dir("test1x","test_no_prev_3",force=True)

    # Create a dataset to copy in
    sv = Supervisor()
    sv.create_calc_dir("test2x","something",df=df_csv,gene_tree=gene_tree)
    shutil.copy(os.path.join(sv.input_dir,"dataframe.csv"),
                os.path.join(sv.output_dir,"dataframe.csv"))
    shutil.copy(os.path.join(sv.input_dir,"gene-tree.newick"),
                os.path.join(sv.output_dir,"gene-tree.newick"))
    sv.finalize()

    # make sure dataframe and tree are properly copied in; other files
    sv.create_calc_dir("test3x","something_else")
    assert os.path.isfile(os.path.join(sv.input_dir,"dataframe.csv"))
    assert os.path.isfile(os.path.join(sv.input_dir,"gene-tree.newick"))
    sv.finalize()

    # Make sure we can load in a model
    assert sv.model is None
    sv.create_calc_dir("test4","something_else_again",model="october")
    assert sv.model == "october"
    sv.finalize()

    # Make sure we can load gene, species, and reconciled trees
    sv.create_calc_dir("test5",
                       calc_type="load it all",
                       df=df_csv,
                       gene_tree=gene_tree,
                       species_tree=species_tree,
                       reconciled_tree=reconciled_tree,
                       model="november")
    assert os.path.isfile(os.path.join(sv.input_dir,"dataframe.csv"))
    assert os.path.isfile(os.path.join(sv.input_dir,"alignment.phy"))
    assert os.path.isfile(os.path.join(sv.input_dir,"gene-tree.newick"))
    assert os.path.isfile(os.path.join(sv.input_dir,"species-tree.newick"))
    assert os.path.isfile(os.path.join(sv.input_dir,"reconciled-tree.newick"))
    assert sv.model == "november"
    sv.finalize()

    new_tree = ete3.Tree("((A,B),(C,D));")
    sv.create_calc_dir("test6",calc_type="yo",gene_tree=new_tree)
    assert os.path.isfile(os.path.join(sv.input_dir,"gene-tree.newick"))
    assert sv.gene_tree == os.path.join(sv.input_dir,"gene-tree.newick")

    os.chdir(current_dir)

def test_Supervisor_check_required(tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    sv.create_calc_dir("test","something")
    pathlib.Path(os.path.join(sv.input_dir,"rocks")).touch()
    pathlib.Path(os.path.join(sv.input_dir,"are")).touch()

    sv.update("this","is")
    sv.update("an","attribute")

    sv.check_required(required_files=["rocks","are"],
                     required_values=["this","an"])

    with pytest.raises(FileNotFoundError):
        sv.check_required(required_files=["rocks","are","not a file"])

    with pytest.raises(ValueError):
        sv.check_required(required_values=["this","an","not a file"])

    os.chdir(current_dir)

def test_Supervisor_copy_output_to_output(tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    sv.create_calc_dir("test1","test")

    for x in ["a","b","c"]:

        f = open(os.path.join(sv.output_dir,f"{x}.newick"),"w")
        f.write("")
        f.close()

        f = open(os.path.join(sv.output_dir,f"{x}.txt"),"w")
        f.write("")
        f.close()

    sv.finalize()

    sv.create_calc_dir("test2","test")
    sv.copy_output_to_output("*.newick")

    old_newick = set([os.path.basename(n) for n in glob.glob(os.path.join("test1","output","*.newick"))])
    new_newick = set([os.path.basename(n) for n in glob.glob(os.path.join(sv.output_dir,"*.newick"))])
    assert old_newick == new_newick
    assert len(old_newick) == 3
    assert len(new_newick) == 3

    old_txt = set([os.path.basename(n) for n in glob.glob(os.path.join("test1","output","*.txt"))])
    new_txt = set([os.path.basename(n) for n in glob.glob(os.path.join(sv.output_dir,"*.txt"))])
    assert old_txt != new_txt
    assert len(old_txt) == 3
    assert len(new_txt) == 0

    # ----

    sv = Supervisor()
    sv.create_calc_dir("test3","test")

    for x in ["a","b","c"]:

        f = open(os.path.join(sv.output_dir,f"{x}.newick"),"w")
        f.write("")
        f.close()

        f = open(os.path.join(sv.output_dir,f"{x}.txt"),"w")
        f.write("")
        f.close()

    sv.finalize()

    sv.create_calc_dir("test4","test")
    sv.copy_output_to_output("*.txt")

    old_newick = set([os.path.basename(n) for n in glob.glob(os.path.join("test3","output","*.newick"))])
    new_newick = set([os.path.basename(n) for n in glob.glob(os.path.join(sv.output_dir,"*.newick"))])
    assert old_newick != new_newick
    assert len(old_newick) == 3
    assert len(new_newick) == 0

    old_txt = set([os.path.basename(n) for n in glob.glob(os.path.join("test3","output","*.txt"))])
    new_txt = set([os.path.basename(n) for n in glob.glob(os.path.join(sv.output_dir,"*.txt"))])
    assert old_txt == new_txt
    assert len(old_txt) == 3
    assert len(new_txt) == 3


    os.chdir(current_dir)

def test_Supervisor_stash(tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    sv.create_calc_dir("stash_0","some_calc")

    # Make some test files and directories
    f = open("stash_a","w")
    f.write("a")
    f.close()

    f = open("stash_b","w")
    f.write("b")
    f.close()

    os.mkdir("stash_c_dir")
    os.mkdir(os.path.join("stash_c_dir","sub_dir"))

    c_file = os.path.join("stash_c_dir","sub_dir","stash_c")
    f = open(c_file,'w')
    f.write('c')
    f.close()

    sv.stash("stash_a")
    assert os.path.isfile(os.path.join(sv.output_dir,"stash_a"))
    sv.stash("stash_b")
    assert os.path.isfile(os.path.join(sv.output_dir,"stash_a"))
    assert os.path.isfile(os.path.join(sv.output_dir,"stash_b"))

    # Test copying in directory
    sv.stash("stash_c_dir")
    assert os.path.isdir(os.path.join(sv.output_dir,"stash_c_dir"))
    assert os.path.isfile(os.path.join(sv.output_dir,"stash_c_dir","sub_dir","stash_c"))

    # Bad file
    with pytest.raises(FileNotFoundError):
        sv.stash("not_a_file")

    # target_dir
    sv.stash("stash_a",target_dir="input")
    sv.stash("stash_b",target_dir="working")
    assert os.path.isfile(os.path.join(sv.input_dir,"stash_a"))
    assert os.path.isfile(os.path.join(sv.working_dir,"stash_b"))

    with pytest.raises(ValueError):
        sv.stash("stash_a",target_dir="not_a_target_dir")

    # Copy in file with target name and make sure it's the right file
    sv.stash("stash_a",target_name="stash_d")
    assert os.path.isfile(os.path.join(sv.output_dir,"stash_d"))
    f = open(os.path.join(sv.output_dir,"stash_d"))
    contents = f.read()
    f.close()
    assert contents.strip() == "a"

    # Send in a subdirectory and make sure the file goes in as expected
    sv.stash("stash_b",target_name="sub_dir_with_file/file.txt")
    assert os.path.isdir(os.path.join(sv.output_dir,"sub_dir_with_file"))
    assert os.path.isfile(os.path.join(sv.output_dir,"sub_dir_with_file","file.txt"))
    f = open(os.path.join(sv.output_dir,"sub_dir_with_file","file.txt"))
    contents = f.read()
    f.close()
    assert contents.strip() == "b"

    # Try to use an existing file as a directory...
    sv.stash("stash_b",target_name="stash_d/file.txt")
    assert os.path.isdir(os.path.join(sv.output_dir,"stash_d"))
    assert os.path.isfile(os.path.join(sv.output_dir,"stash_d","file.txt"))

    os.chdir(current_dir)

def test_Supervisor_write_json():

    # Tested throughout the function. As of right now it's super simple, so
    # not testing specifically
    True

def test_Supervisor_event(tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()

    # --------------------------------------------------------------------------
    # Check for 'running' status.

    assert sv.status == "empty"
    with pytest.raises(RuntimeError):
        sv.event("event")

    sv._run_parameters["calc_status"] = "complete"
    with pytest.raises(RuntimeError):
        sv.event("event")

    sv._run_parameters["calc_status"] = "crashed"
    with pytest.raises(RuntimeError):
        sv.event("event")

    sv = Supervisor()
    sv.create_calc_dir("test0.0","test0.0")
    assert sv.status == "running"
    sv.event("event")

    # --------------------------------------------------------------------------
    # Make sure no events recorded when it starts
    sv = Supervisor()
    sv.create_calc_dir("test0","test0")
    assert len(sv.run_parameters["events"]) == 0

    # --------------------------------------------------------------------------
    # Send in reserved keywords
    with pytest.raises(ValueError):
        sv.event("test",local_directory='5')

    # Send in reserved keywords
    with pytest.raises(ValueError):
        sv.event("test",time=5)

    # --------------------------------------------------------------------------
    # Basic event recording
    sv = Supervisor()
    sv.create_calc_dir("test1","test1")
    sv.event("first_event",keyword="saturn v")
    assert len(sv.run_parameters["events"]) == 1
    evt = sv.run_parameters["events"][0]
    assert evt["description"] == "first_event"
    assert evt["local_directory"] == os.getcwd()
    assert evt["time"] <= time.time()
    assert evt["keyword"] == "saturn v"


    # Change directories and make sure relative directories are stored correctly
    os.chdir("test1")
    sv.event("second_event",keyword="atlas")
    assert len(sv.run_parameters["events"]) == 2
    evt = sv.run_parameters["events"][1]
    assert evt["description"] == "second_event"
    assert evt["local_directory"] == "."
    assert evt["time"] <= time.time()
    assert evt["keyword"] == "atlas"

    # Change directories and see what happens
    os.chdir("working")
    sv.event("third_event",keyword="gemini")
    assert len(sv.run_parameters["events"]) == 3
    evt = sv.run_parameters["events"][2]
    assert evt["description"] == "third_event"
    assert evt["local_directory"] == "working"
    assert evt["time"] <= time.time()
    assert evt["keyword"] == "gemini"

    # --------------------------------------------------------------------------
    # Make sure this records to run_paramters.json as it is run

    sv = Supervisor()
    sv.create_calc_dir("test2","test2")
    sv.event("first_event",keyword="saturn v")
    assert len(sv.run_parameters["events"]) == 1
    evt = sv.run_parameters["events"][0]
    assert evt["description"] == "first_event"
    assert evt["local_directory"] == os.getcwd()
    assert evt["time"] <= time.time()
    assert evt["keyword"] == "saturn v"

    json_file = os.path.join(sv.calc_dir,"run_parameters.json")
    f = open(json_file)
    p = json.load(f)
    f.close()

    assert p["events"][0]["description"] == "first_event"
    assert p["events"][0]["local_directory"] == evt["local_directory"]
    assert p["events"][0]["time"] == evt["time"]
    assert p["events"][0]["keyword"] == evt["keyword"]

    os.chdir(current_dir)

def test_Supervisor_finalize(tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()

    # should ignore finalize with empty
    sv.finalize()
    assert sv.status == "empty"

    sv.create_calc_dir("stupid","another")
    sv.finalize()
    assert sv.status == "complete"
    json_file = os.path.join(sv.calc_dir,"run_parameters.json")
    f = open(json_file)
    p = json.load(f)
    f.close()
    assert p["calc_status"] == "complete"
    assert p["completion_time"] <= time.time()

    sv = Supervisor()
    sv.create_calc_dir("stupid2","another")

    os.chdir("stupid2")
    assert os.getcwd() != sv.starting_dir
    sv.finalize(successful=False)
    assert sv.status == "crashed"
    json_file = os.path.join(sv.calc_dir,"run_parameters.json")
    f = open(json_file)
    p = json.load(f)
    f.close()
    assert p["calc_status"] == "crashed"
    assert p["completion_time"] <= time.time()

    # This tests return-to-starting directory after crash functionality
    assert os.getcwd() == sv.starting_dir

    os.chdir(current_dir)

def test_Supervisor_alignment(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.alignment is None

    sv.create_calc_dir("test0","test0",df=simple_phylo["dataframe.csv"])
    assert sv.alignment == os.path.join(sv.input_dir,"alignment.phy")

    os.chdir(current_dir)

def test_Supervisor_df(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.df is None

    sv.create_calc_dir("test0","test0",df=simple_phylo["dataframe.csv"])
    assert issubclass(type(sv.df),pd.DataFrame)

    os.chdir(current_dir)

def test_Supervisor_gene_tree(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.gene_tree is None

    sv.create_calc_dir("test0","test0",gene_tree=simple_phylo["tree.newick"])
    assert sv.gene_tree == os.path.join(sv.input_dir,"gene-tree.newick")

    os.chdir(current_dir)

def test_Supervisor_species_tree(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.species_tree is None

    # Just testing load, not type of tree, so okay to bring in this tree
    sv.create_calc_dir("test0","test0",species_tree=simple_phylo["tree.newick"])
    assert sv.species_tree == os.path.join(sv.input_dir,"species-tree.newick")

    os.chdir(current_dir)

def test_Supervisor_reconciled_tree(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.reconciled_tree is None

    # Just testing load, not type of tree, so okay to bring in this tree
    sv.create_calc_dir("test0","test0",reconciled_tree=simple_phylo["tree.newick"])
    assert sv.reconciled_tree == os.path.join(sv.input_dir,"reconciled-tree.newick")

    os.chdir(current_dir)

def test_Supervisor_calc_dir(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.calc_dir is None

    sv.create_calc_dir("test0","test0")
    assert sv.calc_dir == os.path.abspath("test0")

    os.chdir(current_dir)

def test_Supervisor_model(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.model is None

    sv.create_calc_dir("test0","test0",model="JTT")
    assert sv.model == "JTT"

    os.chdir(current_dir)

def test_Supervisor_calc_type(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.calc_type is None

    sv.create_calc_dir("test0",calc_type="testy")
    assert sv.calc_type == "testy"

    os.chdir(current_dir)


def test_Supervisor_previous_entries(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.previous_entries is None

    sv.create_calc_dir("test0",calc_type="test0")
    assert sv.previous_entries is None
    sv.finalize()

    sv.create_calc_dir("test1",calc_type="test1")
    assert len(sv.previous_entries) == 1
    assert sv.previous_entries[0]["calc_dir"] == os.path.abspath("test0")

    os.chdir(current_dir)

def test_Supervisor_run_parameters(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert len(sv.run_parameters) == 4
    assert isinstance(sv.run_parameters,dict)

    sv.create_calc_dir("test0",calc_type="testy")
    assert len(sv.run_parameters) > 4

    os.chdir(current_dir)

def test_Supervisor_seed(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert isinstance(sv.seed,int)
    assert sv.seed > 0
    this_seed = sv.seed

    sv.create_calc_dir("test0",calc_type="testy")
    assert this_seed == sv.seed

    os.chdir(current_dir)


def test_Supervisor_status(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.status == "empty"

    sv.create_calc_dir("test0",calc_type="testy")
    assert sv.status == "running"

    os.chdir(current_dir)

def test_Supervisor_update(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()

    with pytest.raises(KeyError):
        sv.run_parameters["rocket"]

    sv.update("rocket","launch")
    assert sv.run_parameters["rocket"] == "launch"

    with pytest.raises(ValueError):
        sv.update((1,2),"stupid")

    with pytest.raises(KeyError):
        sv.run_parameters["another_key"]

    sv.create_calc_dir("test0","test0")
    sv.update("another_key",(1,2))
    assert sv.run_parameters["another_key"] == (1,2)

    sv.finalize()
    with pytest.raises(RuntimeError):
        sv.update("key_after_finalize",False)

    sv._run_parameters["calc_status"] = "running"
    sv.update("okay_because_running",1)

    sv.finalize(successful=False)
    with pytest.raises(RuntimeError):
        sv.update("key_after_finalize",False)

    os.chdir(current_dir)

def test_Supervisor_input_dir(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.input_dir is None

    sv.create_calc_dir("test0","test0")
    assert sv.input_dir == os.path.abspath(os.path.join("test0","input"))

    os.chdir(current_dir)


def test_Supervisor_output_dir(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.output_dir is None

    sv.create_calc_dir("test0","test0")
    assert sv.output_dir == os.path.abspath(os.path.join("test0","output"))

    os.chdir(current_dir)

def test_Supervisor_starting_dir(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.starting_dir == os.path.abspath(tmpdir)

    sv.create_calc_dir("test0","test0")
    assert sv.starting_dir == os.path.abspath(tmpdir)

    os.chdir(current_dir)

def test_Supervisor_working_dir(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.working_dir is None

    sv.create_calc_dir("test0","test0")
    assert sv.working_dir == os.path.abspath(os.path.join("test0","working"))

    os.chdir(current_dir)

def test_Supervisor_tree_class(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    sv = Supervisor()
    assert sv.tree_class is None

    sv.create_calc_dir("test0","ml_tree")
    assert sv.tree_class == "gene"
    sv.finalize()

    sv.create_calc_dir("test2","reconcile_tree")
    assert sv.tree_class == "reconciled"
    sv.finalize()

    sv.create_calc_dir("test3","something")
    assert sv.tree_class is None
    sv.finalize()

    os.chdir(current_dir)
