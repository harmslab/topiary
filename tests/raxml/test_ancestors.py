import pytest

import topiary
from topiary.raxml.ancestors import _make_ancestor_summary_trees
from topiary.raxml.ancestors import _parse_raxml_anc_output
from topiary.raxml.ancestors import generate_ancestors
from topiary.raxml import RAXML_BINARY

import numpy as np

import os, ete3

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test__make_ancestor_summary_trees():

    pass

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test__parse_raxml_anc_output():

    pass

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_generate_ancestors(simple_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs = {"previous_dir":None,
              "df":simple_phylo["dataframe.csv"],
              "model":"JTT",
              "tree_file":simple_phylo["tree.newick"],
              "alt_cutoff":0.25,
              "calc_dir":"ancestors",
              "overwrite":False,
              "supervisor":None,
              "num_threads":1,
              "raxml_binary":RAXML_BINARY}

    generate_ancestors(**kwargs)

    output = os.path.join("ancestors","output")

    expected = ["tree_anc-label.newick",
                "tree_anc-pp.newick",
                "summary-tree.pdf",
                "dataframe.csv"]

    expected.append(os.path.join("ancestors","ancestors.fasta"))
    expected.append(os.path.join("ancestors","ancestor-data.csv"))
    for i in range(6):
        expected.append(os.path.join("ancestors",f"anc{i+1}.pdf"))

    for e in expected:
        assert os.path.isfile(os.path.join(output,e))
