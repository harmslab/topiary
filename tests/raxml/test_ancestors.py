import pytest

import topiary
from topiary.raxml.ancestors import _make_ancestor_summary_trees
from topiary.raxml.ancestors import _parse_raxml_anc_output
from topiary.raxml.ancestors import generate_ancestors
from topiary.raxml import RAXML_BINARY

import ete3

import numpy as np

import os
import copy
import json

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test__make_ancestor_summary_trees():

    pass

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test__parse_raxml_anc_output():

    pass

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_generate_ancestors(tiny_phylo,tmpdir):

    df = tiny_phylo["initial-input/dataframe.csv"]
    gene_tree = tiny_phylo["final-output/gene-tree.newick"]
    species_tree = tiny_phylo["initial-input/species-tree.newick"]
    reconciled_tree = tiny_phylo["final-output/reconciled-tree.newick"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs_template = {"prev_calculation":None,
                       "df":df,
                       "model":"JTT",
                       "gene_tree":gene_tree,
                       "alt_cutoff":0.25,
                       "calc_dir":"ancestors",
                       "overwrite":False,
                       "num_threads":1,
                       "raxml_binary":RAXML_BINARY}

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["calc_dir"] = "test0"

    generate_ancestors(**kwargs)

    output = os.path.join("test0","output")

    expected = ["gene-tree_anc-label.newick",
                "gene-tree_anc-pp.newick",
                "summary-tree.pdf",
                "dataframe.csv"]

    expected.append(os.path.join("gene-tree_ancestors","ancestors.fasta"))
    expected.append(os.path.join("gene-tree_ancestors","ancestor-data.csv"))
    for i in range(6):
        expected.append(os.path.join("gene-tree_ancestors",f"anc{i+1}.pdf"))

    for e in expected:
        assert os.path.isfile(os.path.join(output,e))

    f = open(os.path.join("test0","run_parameters.json"))
    run_params = json.load(f)
    f.close()
    assert run_params["model"] == "JTT"
    assert run_params["alt_cutoff"] == 0.25



    # --------------------------------------------------------------------------
    # Make sure reconciled tree takes precendence over gene tree

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["gene_tree"] = gene_tree
    kwargs["reconciled_tree"] = reconciled_tree
    kwargs["model"] = "LG"
    kwargs["alt_cutoff"] = 0.20
    kwargs["calc_dir"] = "test1"

    generate_ancestors(**kwargs)

    output = os.path.join("test1","output")

    expected = ["reconciled-tree_anc-label.newick",
                "reconciled-tree_anc-pp.newick",
                "summary-tree.pdf",
                "dataframe.csv"]

    expected.append(os.path.join("reconciled-tree_ancestors","ancestors.fasta"))
    expected.append(os.path.join("reconciled-tree_ancestors","ancestor-data.csv"))
    for i in range(6):
        expected.append(os.path.join("reconciled-tree_ancestors",f"anc{i+1}.pdf"))

    for e in expected:
        assert os.path.isfile(os.path.join(output,e))

    f = open(os.path.join("test1","run_parameters.json"))
    run_params = json.load(f)
    f.close()
    assert run_params["model"] == "LG"
    assert run_params["alt_cutoff"] == 0.20

    os.chdir("..")

    # --------------------------------------------------------------------------
    # make sure it dies when neither gene tree nor reconciled passed in

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["calc_dir"] = "test2"
    kwargs["gene_tree"] = None
    kwargs["reconciled_tree"] = None

    with pytest.raises(ValueError):
        generate_ancestors(**kwargs)
