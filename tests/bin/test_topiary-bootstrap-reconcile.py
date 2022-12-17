import pytest

import topiary

import os
import subprocess
import shutil
import json
import glob

import tarfile

@pytest.mark.run_generax
@pytest.mark.run_raxml
def test_main(tiny_phylo,tmpdir):

    prev_bs = tiny_phylo["04_bootstraps_toy"]
    current_bs = tiny_phylo["05_reconcile-bootstraps_toy"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    # Get location of binary
    location = os.path.dirname(os.path.realpath(__file__))
    test_bin = os.path.join(location,"..","..","bin","topiary-bootstrap-reconcile")
    base_cmd = [test_bin]

    # Should run but fail because no arguments
    ret = subprocess.run(base_cmd)
    assert ret.returncode != 0

    print("----> TEST 0")

    os.mkdir("test0")
    os.chdir("test0")
    os.mkdir("existing-run")
    shutil.copytree(prev_bs,os.path.join("existing-run","04_bootstraps"))
    cmd = base_cmd[:]
    cmd.append("existing-run")

    # Fail because no threads
    ret = subprocess.run(cmd)
    assert ret.returncode != 0

    # should fail because more threads than replicates
    cmd.append("10")
    ret = subprocess.run(cmd)
    assert ret.returncode != 0

    # This should run
    cmd[-1] = "2"

    ret = subprocess.run(cmd)
    assert ret.returncode == 0
    assert os.getcwd() == os.path.join(tmpdir,"test0")

    out_base = os.path.join("existing-run","05_reconciled-tree-bootstraps")
    expected_files = ["dataframe.csv",
                      "reconciled-tree.newick",
                      "reconciled-tree_events.newick",
                      "species-tree.newick",
                      "gene-tree.newick",
                      "reconciled-tree_anc-label.newick",
                      "reconciled-tree_supports.newick",
                      "summary-tree.pdf",
                      "gene-tree_supports.newick",
                      "reconciled-tree_anc-pp.newick",
                      "reconciliations.txt"]
    for e in expected_files:
        assert os.path.isfile(os.path.join(out_base,"output",e))

    json_file = os.path.join(out_base,"run_parameters.json")
    f = open(json_file,"r")
    param = json.load(f)
    f.close()
    assert param["calc_type"] == "reconcile_bootstrap"
    assert param["calc_status"] == "complete"
    assert param["model"] == "LG"
    assert param["allow_horizontal_transfer"] == True
    assert param["seed"] == 12345
    assert param["bootstrap_converged"] == False

    os.chdir("..")

    # print("----> TEST 1")
    #
    # # Test restart capability
    #
    # os.mkdir("test1")
    # os.chdir("test1")
    # os.mkdir("existing-run")
    # shutil.copytree(prev_bs,os.path.join("existing-run","04_bootstraps"))
    # shutil.copytree(os.path.join("..","test0","existing-run","05_reconcile-bootstraps"),
    #                 os.path.join("existing-run","05_reconcile-bootstraps"))
    # cmd = base_cmd[:]
    # cmd.append("existing-run")
    # cmd.append("2")
    #
    # # Set status to running rather than complete
    # f = open(os.path.join("existing-run","05_reconcile-bootstraps","run_parameters.json"))
    # run_param = json.load(f)
    # f.close()
    # run_param["calc_status"] = "running"
    # f = open(os.path.join("existing-run","05_reconcile-bootstraps","run_parameters.json"),"w")
    # json.dump(run_param,f)
    # f.close()
    #
    # working_dir = os.path.join("existing-run","05_reconcile-bootstraps","working")
    # os.chdir(working_dir)
    #
    # # open file
    # f = tarfile.open('replicates.tar.gz')
    # f.extractall('replicates')
    # f.close()
    #
    # os.chdir("../../../")
    #
    # # Wipe out some claim files
    # to_nuke = glob.glob(os.path.join(working_dir,"replicates","00003","complete*"))
    # for n in to_nuke:
    #     os.remove(n)
    # to_nuke = glob.glob(os.path.join(working_dir,"replicates","00003","claim*"))
    # for n in to_nuke:
    #     os.remove(n)
    # to_nuke = glob.glob(os.path.join(working_dir,"replicates","00004","complete*"))
    # for n in to_nuke:
    #     os.remove(n)
    # to_nuke = glob.glob(os.path.join(working_dir,"replicates","00004","claim*"))
    # for n in to_nuke:
    #     os.remove(n)
    #
    # # Make sure output bit is gone so we can see if it ran
    # os.remove(os.path.join("existing-run","05_reconcile-bootstraps","output","reconciled-tree_supports.newick"))
    #
    # # Should not work because neither --overwrite nor --restart set
    # ret = subprocess.run(cmd)
    # assert ret.returncode == 1
    #
    # # Should work and go down restart route. No easy way to check to see if
    # # that's what it did other than seeing that it runs.
    # cmd.append("--restart")
    # ret = subprocess.run(cmd)
    # assert ret.returncode == 0
    #
    # assert os.getcwd() == os.path.join(tmpdir,"test1")
    #
    # out_base = os.path.join("existing-run","05_reconcile-bootstraps")
    # expected_files = ["dataframe.csv",
    #                   "reconciled-tree.newick",
    #                   "reconciled-tree_events.newick",
    #                   "species-tree.newick",
    #                   "gene-tree.newick",
    #                   "reconciled-tree_anc-label.newick",
    #                   "reconciled-tree_supports.newick",
    #                   "summary-tree.pdf",
    #                   "gene-tree_supports.newick",
    #                   "reconciled-tree_anc-pp.newick",
    #                   "reconciliations.txt"]
    # for e in expected_files:
    #     assert os.path.isfile(os.path.join(out_base,"output",e))
    #
    # json_file = os.path.join(out_base,"run_parameters.json")
    # f = open(json_file,"r")
    # param = json.load(f)
    # f.close()
    # assert param["calc_type"] == "reconcile_bootstrap"
    # assert param["calc_status"] == "complete"
    # assert param["model"] == "LG"
    # assert param["allow_horizontal_transfer"] == True
    # assert param["seed"] == 12345
    # assert param["bootstrap_converged"] == False
    #
    # os.chdir("..")
    #
    # print("----> TEST 2")
    #
    # # Test overwrite flag
    #
    # os.mkdir("test2")
    # os.chdir("test2")
    # os.mkdir("existing-run")
    # shutil.copytree(prev_bs,os.path.join("existing-run","04_bootstraps"))
    # shutil.copytree(os.path.join("..","test0","existing-run","05_reconcile-bootstraps"),
    #                 os.path.join("existing-run","05_reconcile-bootstraps"))
    # cmd = base_cmd[:]
    # cmd.append("existing-run")
    # cmd.append('2')
    # cmd.append("--overwrite")
    #
    # # Make sure output bit is gone so we can see if it ran
    # os.remove(os.path.join("existing-run","05_reconcile-bootstraps","output","reconciled-tree_supports.newick"))
    #
    # ret = subprocess.run(cmd)
    # assert ret.returncode == 0
    # assert os.getcwd() == os.path.join(tmpdir,"test0")
    #
    # out_base = os.path.join("existing-run","05_reconcile-bootstraps")
    # expected_files = ["dataframe.csv",
    #                   "reconciled-tree.newick",
    #                   "reconciled-tree_events.newick",
    #                   "species-tree.newick",
    #                   "gene-tree.newick",
    #                   "reconciled-tree_anc-label.newick",
    #                   "reconciled-tree_supports.newick",
    #                   "summary-tree.pdf",
    #                   "gene-tree_supports.newick",
    #                   "reconciled-tree_anc-pp.newick",
    #                   "reconciliations.txt"]
    # for e in expected_files:
    #     assert os.path.isfile(os.path.join(out_base,"output",e))
    #
    # json_file = os.path.join(out_base,"run_parameters.json")
    # f = open(json_file,"r")
    # param = json.load(f)
    # f.close()
    # assert param["calc_type"] == "reconcile_bootstrap"
    # assert param["calc_status"] == "complete"
    # assert param["model"] == "LG"
    # assert param["allow_horizontal_transfer"] == True
    # assert param["seed"] == 12345
    # assert param["bootstrap_converged"] == False
    #
    # os.chdir("..")


    os.chdir(current_dir)
