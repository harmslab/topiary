import pytest
import topiary

from topiary.raxml.bootstrap import generate_bootstraps
from topiary.raxml import RAXML_BINARY

import os
import json

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_generate_bootstraps(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs = {"previous_dir":None,
              "df":simple_phylo["dataframe.csv"],
              "model":"JTT",
              "gene_tree":simple_phylo["tree_ml.newick"],
              "calc_dir":"bootstraps",
              "overwrite":False,
              "num_bootstraps":10,
              "supervisor":None,
              "num_threads":1,
              "raxml_binary":RAXML_BINARY}

    generate_bootstraps(**kwargs)

    expected_files = ["summary-tree.pdf",
                      "gene-tree_supports.newick",
                      "dataframe.csv"]
    for e in expected_files:
        assert os.path.isfile(os.path.join("bootstraps","output",e))

    json_file = os.path.join("bootstraps","run_parameters.json")
    assert os.path.isfile(json_file)
    f = open(json_file,"r")
    param = json.load(f)
    f.close()

    assert param["calc_type"] == "ml_bootstrap"
    assert param["model"] == "JTT"

    os.chdir(current_dir)
