import pytest
import topiary

from topiary.generax.reconcile import reconcile
from topiary.generax import GENERAX_BINARY
from topiary.raxml import RAXML_BINARY
from topiary._private import Supervisor

import os
import copy
import json
import shutil

@pytest.mark.run_generax
def test_reconcile(tiny_phylo,tmpdir):

    df_csv = tiny_phylo["initial-input/dataframe.csv"]
    df = topiary.read_dataframe(df_csv)
    gene_tree = tiny_phylo["final-output/gene-tree.newick"]
    gene_tree_wrong = tiny_phylo["final-output/gene-tree_wrong.newick"]
    species_tree = tiny_phylo["initial-input/species-tree.newick"]
    reconciled_tree = tiny_phylo["final-output/reconciled-tree.newick"]
    bootstrap_directory = tiny_phylo["04_bootstraps_toy/output/bootstrap_replicates"]
    f = open(tiny_phylo["model.txt"],"r")
    model = f.read().strip()
    f.close()

    prev_ml = tiny_phylo["01_gene-tree"]
    prev_bs = tiny_phylo["04_bootstraps_toy"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)


    kwargs_template = {"prev_calculation":None,
                       "df":None,
                       "model":None,
                       "gene_tree":None,
                       "species_tree":None,
                       "reconciled_tree":None,
                       "allow_horizontal_transfer":True,
                       "seed":None,
                       "bootstrap":False,
                       "calc_dir":"reconcile",
                       "overwrite":False,
                       "num_threads":1,
                       "generax_binary":GENERAX_BINARY,
                       "raxml_binary":RAXML_BINARY}

    # Test bad generax binary
    os.mkdir("test0")
    os.chdir("test0")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["generax_binary"] = "not_a_binary"
    with pytest.raises(ValueError):
        reconcile(**kwargs)
    os.chdir("..")

    # Test bad generax binary
    os.mkdir("test0.0")
    os.chdir("test0.0")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["raxml_binary"] = "not_a_binary"
    with pytest.raises(ValueError):
        reconcile(**kwargs)
    os.chdir("..")

    # Test MPI threads check
    os.mkdir("test1")
    os.chdir("test1")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["num_threads"] = 99999999999999
    with pytest.raises(RuntimeError):
        reconcile(**kwargs)
    os.chdir("..")


    # Make sure code looks for model
    os.mkdir("test2")
    os.chdir("test2")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["gene_tree"] = gene_tree
    kwargs["df"] = df_csv
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
    kwargs["df"] = df_csv
    kwargs["model"] = model
    with pytest.raises(FileNotFoundError):
        reconcile(**kwargs)
    os.chdir("..")

    # -------------------------------------------------------------------------
    # No previous calc; no bootstrap

    # Non boostrap run
    os.mkdir("test5.0")
    os.chdir("test5.0")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["gene_tree"] = gene_tree
    kwargs["df"] = df_csv
    kwargs["model"] = model
    kwargs["species_tree"] = species_tree
    kwargs["allow_horizontal_transfer"] = None
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
    assert param["model"] == model
    assert param["allow_horizontal_transfer"] == True

    os.chdir("..")

    # Re-run, set allow_horizontal_transfer False
    os.mkdir("test5.1")
    os.chdir("test5.1")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["gene_tree"] = gene_tree
    kwargs["df"] = df_csv
    kwargs["model"] = model
    kwargs["species_tree"] = species_tree
    kwargs["allow_horizontal_transfer"] = False
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
    assert param["model"] == model
    assert param["allow_horizontal_transfer"] == False

    os.chdir("..")

    # -------------------------------------------------------------------------
    # Previous calc, as directory

    # Non boostrap run
    os.mkdir("test6.0")
    os.chdir("test6.0")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["prev_calculation"] = prev_ml
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
    assert param["model"] == model
    assert param["allow_horizontal_transfer"] == True

    os.chdir("..")

    # Show we can overwrite what's in prev_dir with model and allow_horizontal_transfer

    # Non boostrap run
    os.mkdir("test6.1")
    os.chdir("test6.1")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["prev_calculation"] = prev_ml
    kwargs["model"] = "JTT"
    kwargs["allow_horizontal_transfer"] = False
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
    assert param["allow_horizontal_transfer"] == False

    os.chdir("..")


    # -------------------------------------------------------------------------
    # Previous calc, as supervisor

    # Non boostrap run
    os.mkdir("test7.0")
    os.chdir("test7.0")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["prev_calculation"] = Supervisor(prev_ml)
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
    assert param["model"] == model
    assert param["allow_horizontal_transfer"] == True

    os.chdir("..")

    # Show we can overwrite what's in prev_dir with model and allow_horizontal_transfer

    # Non boostrap run
    os.mkdir("test7.1")
    os.chdir("test7.1")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["prev_calculation"] = Supervisor(prev_ml)
    kwargs["model"] = "JTT"
    kwargs["allow_horizontal_transfer"] = False
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
    assert param["allow_horizontal_transfer"] == False

    os.chdir("..")

    # -------------------------------------------------------------------------
    # Previous calc, bootstrap

    os.mkdir("test8.0")
    os.chdir("test8.0")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["prev_calculation"] = Supervisor(prev_bs)
    kwargs["bootstrap"] = True
    reconcile(**kwargs)

    out_dir = os.path.join("reconcile","output")
    expected_files = ["dataframe.csv",
                      "reconciled-tree.newick",
                      "species-tree.newick",
                      "reconciled-tree_events.newick",
                      "reconciled-tree_supports.newick",
                      "summary-tree.pdf"]
    for e in expected_files:
        assert os.path.isfile(os.path.join(out_dir,e))

    json_file = os.path.join("reconcile","run_parameters.json")
    f = open(json_file,"r")
    param = json.load(f)
    f.close()
    assert param["calc_type"] == "reconcile_bootstrap"
    assert param["model"] == model
    assert param["allow_horizontal_transfer"] == True
    # This should be set
    param["bootstrap_converged"]

    os.chdir("..")

    # Try to bootstrap without a previous bootstrap
    os.mkdir("test8.1")
    os.chdir("test8.1")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["prev_calculation"] = Supervisor(prev_ml)
    kwargs["bootstrap"] = True
    with pytest.raises(ValueError):
        reconcile(**kwargs)

    os.chdir("..")

    # Try to bootstrap but with mangled replicate directory
    os.mkdir("test8.2")
    os.chdir("test8.2")
    shutil.copytree(prev_bs,"mangled-bs")
    shutil.rmtree(os.path.join("mangled-bs","output","bootstrap_replicates"))
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["prev_calculation"] = "mangled-bs"
    kwargs["bootstrap"] = True
    with pytest.raises(FileNotFoundError):
        reconcile(**kwargs)

    os.chdir("..")

    # Try to bootstrap but with mangled reconciled-tree file
    os.mkdir("test8.3")
    os.chdir("test8.3")
    shutil.copytree(prev_bs,"mangled-bs")
    os.remove(os.path.join("mangled-bs","output","reconciled-tree.newick"))
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["prev_calculation"] = "mangled-bs"
    kwargs["bootstrap"] = True
    with pytest.raises(FileNotFoundError):
        reconcile(**kwargs)

    os.chdir("..")

    os.chdir(current_dir)
