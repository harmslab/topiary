import pytest
import topiary

from topiary.raxml.bootstrap import generate_bootstraps
from topiary.raxml import RAXML_BINARY

import os
import json

@pytest.mark.run_raxml
def test_generate_bootstraps(tiny_phylo,tmpdir):

    df = tiny_phylo["initial-input/dataframe.csv"]
    gene_tree = tiny_phylo["final-output/gene-tree.newick"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs = {"prev_calculation":None,
              "df":df,
              "model":"JTT",
              "gene_tree":gene_tree,
              "calc_dir":"bootstraps",
              "overwrite":False,
              "num_bootstraps":10,
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
