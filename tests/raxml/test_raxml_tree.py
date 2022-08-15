import pytest

import topiary
from topiary.raxml.tree import generate_ml_tree
from topiary.raxml._raxml import RAXML_BINARY
from topiary._private import Supervisor

import os
import json

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_generate_ml_tree(tiny_phylo,tmpdir):

    df = tiny_phylo["initial-input/dataframe.csv"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs = {"previous_dir":None,
              "df":df,
              "model":"JTT",
              "calc_dir":"test0",
              "overwrite":False,
              "bootstrap":False,
              "supervisor":None,
              "num_threads":1,
              "raxml_binary":RAXML_BINARY}

    generate_ml_tree(**kwargs)

    expected_files = ["dataframe.csv",
                      "summary-tree.pdf",
                      "gene-tree.newick"]
    for e in expected_files:
        assert os.path.isfile(os.path.join("test0","output",e))

    json_file = os.path.join("test0","run_parameters.json")
    assert os.path.isfile(json_file)
    f = open(json_file,"r")
    param = json.load(f)
    f.close()

    assert param["calc_type"] == "ml_tree"
    assert param["model"] == "JTT"

    supervisor = Supervisor()

    kwargs = {"previous_dir":None,
              "df":df,
              "model":"LG",
              "calc_dir":"test1",
              "overwrite":False,
              "bootstrap":False,
              "supervisor":supervisor,
              "num_threads":1,
              "raxml_binary":RAXML_BINARY}

    generate_ml_tree(**kwargs)

    assert supervisor.starting_dir == os.path.abspath(os.getcwd())
    assert supervisor.calc_dir == os.path.abspath("test1")
    expected_files = ["dataframe.csv","summary-tree.pdf","gene-tree.newick"]
    for e in expected_files:
        assert os.path.isfile(os.path.join("test1","output",e))

    json_file = os.path.join("test1","run_parameters.json")
    assert os.path.isfile(json_file)
    f = open(json_file,"r")
    param = json.load(f)
    f.close()

    assert param["calc_type"] == "ml_tree"
    assert param["model"] == "LG"


    supervisor = Supervisor()

    kwargs = {"previous_dir":None,
              "df":df,
              "model":"LG",
              "calc_dir":"test2",
              "overwrite":False,
              "bootstrap":True,
              "supervisor":supervisor,
              "num_threads":1,
              "raxml_binary":RAXML_BINARY}

    generate_ml_tree(**kwargs)

    assert supervisor.starting_dir == os.path.abspath(os.getcwd())
    assert supervisor.calc_dir == os.path.abspath("test2")
    expected_files = ["dataframe.csv",
                      "summary-tree.pdf",
                      "gene-tree.newick",
                      "gene-tree_supports.newick"]
    for e in expected_files:
        assert os.path.isfile(os.path.join("test2","output",e))

    assert os.path.isdir(os.path.join("test2","output","bootstrap_replicates"))
    assert os.path.isfile(os.path.join("test2","output","bootstrap_replicates","bootstraps.newick"))

    json_file = os.path.join("test2","run_parameters.json")
    assert os.path.isfile(json_file)
    f = open(json_file,"r")
    param = json.load(f)
    f.close()

    assert param["calc_type"] == "ml_bootstrap"
    assert param["model"] == "LG"

    os.chdir(current_dir)
