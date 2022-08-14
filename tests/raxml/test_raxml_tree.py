import pytest

import topiary
from topiary.raxml.tree import generate_ml_tree
from topiary.raxml._raxml import RAXML_BINARY
from topiary._private import Supervisor

import os
import json

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_generate_ml_tree(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs = {"previous_dir":None,
              "df":simple_phylo["dataframe.csv"],
              "model":"JTT",
              "calc_dir":"ml_tree_0",
              "overwrite":False,
              "bootstrap":False,
              "supervisor":None,
              "num_threads":1,
              "raxml_binary":RAXML_BINARY}

    generate_ml_tree(**kwargs)

    expected_files = ["dataframe.csv",
                      "summary-tree.pdf",
                      "tree.newick"]
    for e in expected_files:
        assert os.path.isfile(os.path.join("ml_tree_0","output",e))

    json_file = os.path.join("ml_tree_0","run_parameters.json")
    assert os.path.isfile(json_file)
    f = open(json_file,"r")
    param = json.load(f)
    f.close()

    assert param["calc_type"] == "ml_tree"
    assert param["model"] == "JTT"

    supervisor = Supervisor()

    kwargs = {"previous_dir":None,
              "df":simple_phylo["dataframe.csv"],
              "model":"LG",
              "calc_dir":"ml_tree",
              "overwrite":False,
              "bootstrap":False,
              "supervisor":supervisor,
              "num_threads":1,
              "raxml_binary":RAXML_BINARY}

    generate_ml_tree(**kwargs)

    assert supervisor.starting_dir == os.path.abspath(os.getcwd())
    assert supervisor.calc_dir == os.path.abspath("ml_tree")
    expected_files = ["dataframe.csv","summary-tree.pdf","tree.newick"]
    for e in expected_files:
        assert os.path.isfile(os.path.join("ml_tree","output",e))

    json_file = os.path.join("ml_tree","run_parameters.json")
    assert os.path.isfile(json_file)
    f = open(json_file,"r")
    param = json.load(f)
    f.close()

    assert param["calc_type"] == "ml_tree"
    assert param["model"] == "LG"

    os.chdir(current_dir)
