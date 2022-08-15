import pytest
import topiary

from topiary.generax.reconcile import reconcile
from topiary.generax import GENERAX_BINARY

import os
import copy
import json

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_reconcile(generax_data,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    input_dir = os.path.abspath(os.path.join(generax_data["toy-input"],"toy-bootstrap","output"))
    df = topiary.read_dataframe(os.path.join(input_dir,"dataframe.csv"))
    model = "JTT"
    gene_tree = os.path.join(os.path.join(generax_data["toy-input"],"toy-ml","output","tree_wrong.newick"))
    species_tree = os.path.join(input_dir,"species_tree.newick")
    bootstrap_directory = os.path.join(input_dir,"bootstrap_replicates")

    kwargs_template = {"previous_dir":None,
                       "df":None,
                       "model":None,
                       "gene_tree":None,
                       "species_tree":None,
                       "allow_horizontal_transfer":True,
                       "bootstrap":False,
                       "calc_dir":"reconcile",
                       "overwrite":False,
                       "supervisor":None,
                       "num_threads":1,
                       "generax_binary":GENERAX_BINARY}

    # Test bad generax binary
    os.mkdir("test0")
    os.chdir("test0")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["generax_binary"] = "not_a_binary"
    with pytest.raises(ValueError):
        reconcile(**kwargs)
    os.chdir("..")

    # Test MPI threads check
    os.mkdir("test1")
    os.chdir("test1")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["num_threads"] = 99999999999999
    with pytest.raises(ValueError):
        reconcile(**kwargs)
    os.chdir("..")

    # Make sure code looks for model
    os.mkdir("test2")
    os.chdir("test2")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["gene_tree"] = gene_tree
    kwargs["df"] = df
    #kwargs["model"] = model
    with pytest.raises(ValueError):
        reconcile(**kwargs)
    os.chdir("..")

    # Make sure code looks for df
    os.mkdir("test3")
    os.chdir("test3")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["gene_tree"] = gene_tree
    #kwargs["df"] = df
    kwargs["model"] = model
    with pytest.raises(FileNotFoundError):
        reconcile(**kwargs)
    os.chdir("..")

    # Make sure code looks for gene_tree
    os.mkdir("test4")
    os.chdir("test4")
    kwargs = copy.deepcopy(kwargs_template)
    #kwargs["gene_tree"] = gene_tree
    kwargs["df"] = df
    kwargs["model"] = model
    with pytest.raises(FileNotFoundError):
        reconcile(**kwargs)
    os.chdir("..")

    # Non boostrap run
    os.mkdir("test5")
    os.chdir("test5")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["gene_tree"] = gene_tree
    kwargs["df"] = df
    kwargs["model"] = model
    kwargs["species_tree"] = species_tree
    reconcile(**kwargs)

    out_dir = os.path.join("reconcile","output")
    expected_files = ["dataframe.csv",
                      "reconciled-tree.newick",
                      "species-tree.newick",
                      "reconciled-tree_events.newick",
                      "summary-tree.pdf"]
    for e in expected_files:
        assert os.path.isfile(os.path.join(out_dir,e))

    json_file = os.path.join("reconcile","run_parameters.json")
    f = open(json_file,"r")
    param = json.load(f)
    f.close()
    assert param["calc_type"] == "reconcile"
    assert param["model"] == "JTT"
    assert param["allow_horizontal_transfer"] == True

    os.chdir("..")

    # Bootstrap run should fail without previous ml run
    os.mkdir("test6")
    os.chdir("test6")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["gene_tree"] = gene_tree
    kwargs["df"] = df
    kwargs["model"] = model
    kwargs["species_tree"] = species_tree
    kwargs["bootstrap"] = True
    with pytest.raises(ValueError):
        reconcile(**kwargs)

    # XX IMPROVE TEST
    # To test this properly I need example outputs. Going to need to get whole
    # pipeline cleaned up to build example.
    # os.mkdir("test7")
    # os.chdir("test7")
    # kwargs = copy.deepcopy(kwargs_template)
    # kwargs["gene_tree"] = gene_tree
    # kwargs["df"] = df
    # kwargs["model"] = "LG"
    # kwargs["species_tree"] = species_tree
    # kwargs["bootstrap"] = True
    # kwargs["allow_horizontal_transfer"] = False
    #
    # reconcile(**kwargs)
    #
    # out_dir = os.path.join("reconcile","output")
    # expected_files = ["dataframe.csv",
    #                   "reconciled-tree.newick",
    #                   "species-tree.newick",
    #                   "reconciled-tree_events.newick",
    #                   "summary-tree.pdf"]
    # for e in expected_files:
    #     assert os.path.isfile(os.path.join(out_dir,e))
    #
    # json_file = os.path.join("reconcile","run_parameters.json")
    # f = open(json_file,"r")
    # param = json.load(f)
    # f.close()
    # assert param["calc_type"] == "reconcile_bootstrap"
    # assert param["model"] == "LG"
    # assert param["allow_horizontal_transfer"] == False
    # param["bootstrap_converged"]


    os.chdir(current_dir)
