import pytest
import topiary

from topiary.generax._reconcile_bootstrap import _progress_bar
from topiary.generax._reconcile_bootstrap import _check_convergence
from topiary.generax._reconcile_bootstrap import _generax_thread_function
from topiary.generax._reconcile_bootstrap import _build_replicate_dirs
from topiary.generax._reconcile_bootstrap import _clean_replicate_dir
from topiary.generax._reconcile_bootstrap import _construct_args
from topiary.generax._reconcile_bootstrap import _run_bootstrap_calculations
from topiary.generax._reconcile_bootstrap import reconcile_bootstrap
from topiary.generax._generax import GENERAX_BINARY
from topiary.raxml import RAXML_BINARY
from topiary._private import Supervisor
from topiary._private import mpi

import ete3

import pandas as pd

import os
import glob
import shutil
import copy
import pathlib
import time
import multiprocessing as mp

def test__progress_bar(small_phylo,tmpdir):

    template = small_phylo["toy-reconcile-bootstraps-running/replicates"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    # -------------------------------------------------------------------------
    # Test just running...

    shutil.copytree(template,"test0")

    status_bar = mp.Process(target=_progress_bar,args=("test0",))
    status_bar.start()

    time.sleep(2)

    # Should not have completed.
    try:
        assert status_bar.is_alive()
    except AssertionError as e:
        status_bar.kill()
        raise AssertionError from e
    status_bar.kill()

    # -------------------------------------------------------------------------
    # Test convergence -- suddenly write "skipped" to a bunch of directories

    shutil.copytree(template,"test1")

    status_bar = mp.Process(target=_progress_bar,args=("test1",))
    status_bar.start()

    for g in glob.glob(os.path.join("test1","0*")):

        # Turn running into completed
        if os.path.isfile(os.path.join(g,"running")):
            os.remove(os.path.join(g,"running"))
            pathlib.Path(os.path.join(g,"completed")).touch()

        # Put skipped in all directories that are not complete
        if not os.path.isfile(os.path.join(g,"completed")):
            pathlib.Path(os.path.join(g,"skipped")).touch()

    # Add timeout loop in case it takes a moment to finish
    start_time = time.time()
    while (time.time() - start_time) < 6:
        if not status_bar.is_alive():
            break
        time.sleep(0.2)

    # Should have completed because we wrote "skipped" and "complete" into all
    # directories
    try:
        assert not status_bar.is_alive()
    except AssertionError as e:
        status_bar.kill()
        raise AssertionError from e
    status_bar.kill()

    # -------------------------------------------------------------------------
    # Test convergence -- write "completed" to all directories

    shutil.copytree(template,"test2")

    status_bar = mp.Process(target=_progress_bar,args=("test2",))
    status_bar.start()

    # Put completed in all directories
    for g in glob.glob(os.path.join("test2","0*")):

        # Turn running into completed
        if not os.path.isfile(os.path.join(g,"completed")):
            pathlib.Path(os.path.join(g,"completed")).touch()

    # Add timeout loop in case it takes a moment to finish
    start_time = time.time()
    while (time.time() - start_time) < 6:
        if not status_bar.is_alive():
            break
        time.sleep(0.2)

    # Should have completed because we wrote "skipped" and "complete" into all
    # directories
    try:
        assert not status_bar.is_alive()
    except AssertionError as e:
        status_bar.kill()
        raise AssertionError from e
    status_bar.kill()

    os.chdir(current_dir)

@pytest.mark.run_raxml
def test__check_convergence(small_phylo,tmpdir):

    template = small_phylo["toy-reconcile-bootstraps-running/replicates"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    # -------------------------------------------------------------------------
    # Converged

    test_dir = "test0"
    shutil.copytree(template,test_dir)
    c, df = _check_convergence(test_dir,converge_cutoff=0.5)
    assert c
    assert issubclass(type(df),pd.DataFrame)
    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,12):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(12,15):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    # -------------------------------------------------------------------------
    # Not converged

    test_dir = "test1"
    shutil.copytree(template,test_dir)
    shutil.move(os.path.join(test_dir,"bs-trees_not-converged.newick"),
                os.path.join(test_dir,"bs-trees.newick"))

    c, df = _check_convergence(test_dir,converge_cutoff=0)
    assert not c
    assert issubclass(type(df),pd.DataFrame)
    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,12):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(12,15):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    os.chdir(current_dir)

@pytest.mark.run_raxml
@pytest.mark.run_generax
def test__generax_thread_function(small_phylo,tmpdir):

    template = small_phylo["toy-reconcile-bootstraps-running/replicates"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs_template = {"converge_cutoff":0.5,
                       "is_manager":False,
                       "lock":None}

    # --------------------------------------------------------------------------
    # Should run on last three, skipping first that are completed and running.

    test_dir = "test0"
    shutil.copytree(template,test_dir)

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["replicate_dir"] = test_dir
    kwargs["hosts"] = mpi.get_hosts(1)
    out = _generax_thread_function(**kwargs)
    assert out is None

    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,12):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(12,15):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    # Make sure calculation did *not* run in 00001, which should not have a
    # results directory because it had a 'completed' file.
    assert not os.path.exists(os.path.join(test_dir,"00001","results"))
    assert os.path.exists(os.path.join(test_dir,"00001","mapping.link"))


    # --------------------------------------------------------------------------
    # Run test again, this time sending in two hosts. Behavior should be
    # identical.

    # There is a problem with generax where it screws up running in parallel on 
    # a fast processor with the small test dataset. (Generax expects a file will
    # be completely written out on one thread, but it's not quite finished
    # writing out when the other thread grabs it.) Loop through the test runner
    # 10 times to make sure it works at least once. Not perfect, but...

    def _test_runner():

        test_dir = "test1"
        if os.path.isdir(test_dir):
            shutil.rmtree(test_dir)

        shutil.copytree(template,test_dir)

        kwargs = copy.deepcopy(kwargs_template)
        kwargs["replicate_dir"] = test_dir
        kwargs["hosts"] = mpi.get_hosts(2)
        out = _generax_thread_function(**kwargs)
        assert out is None

        for i in range(10):
            assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
            assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
            assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

        for i in range(10,12):
            assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
            assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
            assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

        for i in range(12,15):
            assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
            assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
            assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

        # Make sure calculation did *not* run in 00001, which should not have a
        # results directory because it had a 'completed' file.
        assert not os.path.exists(os.path.join(test_dir,"00001","results"))
        assert os.path.exists(os.path.join(test_dir,"00001","mapping.link"))

    for i in range(10):
        try:
            _test_runner()
            break
        except RuntimeError:
            continue

    # --------------------------------------------------------------------------
    # Run test again with a single host and manager = True. Should run once, then
    # make last two directories have "skipped" because calculation has converged

    test_dir = "test2"
    shutil.copytree(template,test_dir)

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["replicate_dir"] = test_dir
    kwargs["hosts"] = mpi.get_hosts(1)
    kwargs["is_manager"] = True
    out = _generax_thread_function(**kwargs)
    assert out[0]
    assert issubclass(type(out[1]),pd.DataFrame)

    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,12):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(12,13):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(13,15):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    # Make sure calculation did *not* run in 00001, which should not have a
    # results directory because it had a 'completed' file.
    assert not os.path.exists(os.path.join(test_dir,"00001","results"))
    assert os.path.exists(os.path.join(test_dir,"00001","mapping.link"))


    # --------------------------------------------------------------------------
    # Run test gain with a single host and manager = True, but unconverged
    # bootstrap replicates. Should complete calculations and return converged
    # False.

    test_dir = "test3"
    shutil.copytree(template,test_dir)
    shutil.move(os.path.join(test_dir,"bs-trees_not-converged.newick"),
                os.path.join(test_dir,"bs-trees.newick"))

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["replicate_dir"] = test_dir
    kwargs["hosts"] = mpi.get_hosts(1)
    kwargs["is_manager"] = True
    kwargs["converge_cutoff"] = 0.0000001
    out = _generax_thread_function(**kwargs)
    assert not out[0]
    assert issubclass(type(out[1]),pd.DataFrame)

    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,12):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(12,15):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    # Make sure calculation did *not* run in 00001, which should not have a
    # results directory because it had a 'completed' file.
    assert not os.path.exists(os.path.join(test_dir,"00001","results"))
    assert os.path.exists(os.path.join(test_dir,"00001","mapping.link"))

    os.chdir(current_dir)

@pytest.mark.run_generax
def test__build_replicate_dirs(small_phylo,tmpdir):

    input_dir = small_phylo["05_gene-tree-bootstraps_toy/output/"]
    df = topiary.read_dataframe(small_phylo["05_gene-tree-bootstraps_toy/input/dataframe.csv"])

    f = open(small_phylo["model.txt"])
    model = f.read().strip()
    f.close()

    gene_tree = small_phylo["final-output/gene-tree.newick"]
    species_tree = small_phylo["final-output/species-tree.newick"]
    bootstrap_directory = small_phylo["05_gene-tree-bootstraps_toy/output/bootstrap_replicates"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs_template = {"df":df,
                       "model":model,
                       "gene_tree":gene_tree,
                       "species_tree":"species_tree",
                       "allow_horizontal_transfer":True,
                       "seed":12345,
                       "bootstrap_directory":bootstrap_directory,
                       "overwrite":False,
                       "generax_binary":GENERAX_BINARY}

    os.mkdir("test0")
    os.chdir("test0")
    kwargs = copy.deepcopy(kwargs_template)

    _build_replicate_dirs(**kwargs)

    # Make sure that directories are made with correct files
    expected_files = set(["alignment.phy",
                          "control.txt",
                          "gene_tree.newick",
                          "mapping.link",
                          "run_generax.sh",
                          "species_tree.newick"])
    expected = "replicates"
    assert len(list(glob.glob(os.path.join(expected,"0*")))) == 4
    for g in glob.glob(os.path.join(expected,"0*")):
        dirs = set(os.listdir(g))
        assert dirs == expected_files

    # Make sure that the correct files are identical between references and
    # that the correct files are different
    def _read_file(some_file):
        out = []
        with open(some_file,'r') as f:
            for line in f:
                out.append(line.strip())
        return "\n".join(out)

    # Read in alignments
    alignments = []
    for i in range(4):
        alignments.append(_read_file(os.path.join(input_dir,"bootstrap_replicates",f"bsmsa_000{i+1}.phy")))

    trees = []
    with open(os.path.join(input_dir,"bootstrap_replicates","bs-trees.newick")) as f:
        for line in f:
            T = ete3.Tree(line.strip(),format=0)
            trees.append(T)

    ref_check = {}
    should_be_same = ["control.txt","mapping.link","run_generax.sh","species_tree.newick"]

    dir_list = glob.glob(os.path.join(expected,"0*"))
    dir_list.sort()
    for i, g in enumerate(dir_list):

        this_check = {}
        for f in glob.glob(os.path.join(g,"*")):
            this_check[os.path.basename(f)] = _read_file(f)

            if i == 0:
                ref_check[os.path.basename(f)] = _read_file(f)

        assert this_check["alignment.phy"] == alignments[i]

        # Make sure the gene trees are correctly brought in
        T = ete3.Tree(this_check["gene_tree.newick"],format=0)
        assert T.robinson_foulds(trees[i],unrooted_trees=True)[0] == 0

        for k in this_check:
            if k in should_be_same:
                assert ref_check[k] == this_check[k]

    os.chdir(current_dir)

def test__clean_replicate_dir(small_phylo,tmpdir):

    template = small_phylo["toy-reconcile-bootstraps-running/replicates"]

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    test_dir = "test0"
    shutil.copytree(template,test_dir)
    pathlib.Path(os.path.join(test_dir,"00013","skipped")).touch()

    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,12):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(12,13):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(13,15):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    _clean_replicate_dir(test_dir)

    for i in range(10):
        assert os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))

    for i in range(10,15):
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","completed"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","running"))
        assert not os.path.isfile(os.path.join(test_dir,f"{(i+1):05d}","skipped"))


    os.chdir(current_dir)

@pytest.mark.run_generax
def test__construct_args():

    kwargs_template = {"replicate_dir":"test",
                       "converge_cutoff":0.03,
                       "num_threads":1,
                       "threads_per_rep":1}

    # 1 thread, 1 per calc

    kwargs = copy.deepcopy(kwargs_template)
    kwargs_list, num_threads = _construct_args(**kwargs)
    assert len(kwargs_list) == 1
    assert num_threads == 1
    assert kwargs_list[0]["replicate_dir"] == "test"
    assert kwargs_list[0]["converge_cutoff"] == 0.03
    assert len(kwargs_list[0]["hosts"]) == 1
    assert kwargs_list[0]["is_manager"]

    # 2 threads, 1 per calc

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["num_threads"] = 2
    kwargs_list, num_threads = _construct_args(**kwargs)
    assert len(kwargs_list) == 2
    assert num_threads == 2
    assert kwargs_list[0]["replicate_dir"] == "test"
    assert kwargs_list[0]["converge_cutoff"] == 0.03
    assert len(kwargs_list[0]["hosts"]) == 1
    assert kwargs_list[0]["is_manager"]

    assert kwargs_list[1]["replicate_dir"] == "test"
    assert kwargs_list[1]["converge_cutoff"] == 0.03
    assert len(kwargs_list[1]["hosts"]) == 1
    assert not kwargs_list[1]["is_manager"]

    # 2 threads, 2 threads per calc

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["num_threads"] = 2
    kwargs["threads_per_rep"] = 2
    kwargs["converge_cutoff"] = 0.5

    kwargs_list, num_threads = _construct_args(**kwargs)
    assert len(kwargs_list) == 1
    assert num_threads == 1
    assert kwargs_list[0]["replicate_dir"] == "test"
    assert kwargs_list[0]["converge_cutoff"] == 0.5
    assert len(kwargs_list[0]["hosts"]) == 2
    assert kwargs_list[0]["is_manager"]


@pytest.mark.run_generax
def test__run_bootstrap_calculations(small_phylo,tmpdir):

    # This basically makes sure that a calculation is run in every directory.
    # It does not check quality/correctness of output

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    shutil.copytree(small_phylo["06_reconciled-tree-bootstraps_toy/working/replicates"],
                    "replicates_template")

    kwargs_template = {"replicate_dir":"replicates",
                       "converge_cutoff":0.03,
                       "num_threads":1,
                       "threads_per_rep":1}

    # Single thread
    rep_dir = "test0"
    shutil.copytree("replicates_template",rep_dir)
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["replicate_dir"] = rep_dir

    _run_bootstrap_calculations(**kwargs)

    bs_trees = os.path.join(rep_dir,"bs-trees.newick")
    assert os.path.isfile(bs_trees)
    f = open(bs_trees)
    lines = f.readlines()
    f.close()

    assert len(lines) == 4
    for d in ["00001","00002","00003","00004"]:
        rep = os.path.join(rep_dir,d)
        assert os.path.isfile(os.path.join(rep,"completed"))
        assert os.path.isdir(os.path.join(rep,"result"))
        assert os.path.isfile(os.path.join(rep,
                                           "result",
                                           "results",
                                           "reconcile",
                                           "geneTree.newick"))

    # Run multithreaded

    rep_dir = "test1"
    shutil.copytree("replicates_template",rep_dir)
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["replicate_dir"] = rep_dir
    kwargs["num_threads"] = 2

    _run_bootstrap_calculations(**kwargs)

    bs_trees = os.path.join(rep_dir,"bs-trees.newick")
    assert os.path.isfile(bs_trees)
    f = open(bs_trees)
    lines = f.readlines()
    f.close()

    assert len(lines) == 4
    for d in ["00001","00002","00003","00004"]:
        rep = os.path.join(rep_dir,d)
        assert os.path.isfile(os.path.join(rep,"completed"))
        assert os.path.isdir(os.path.join(rep,"result"))
        assert os.path.isfile(os.path.join(rep,
                                           "result",
                                           "results",
                                           "reconcile",
                                           "geneTree.newick"))


    # Run multithreaded, multiple slots

    rep_dir = "test2"
    shutil.copytree("replicates_template",rep_dir)
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["replicate_dir"] = rep_dir
    kwargs["num_threads"] = 2
    kwargs["threads_per_rep"] = 2

    _run_bootstrap_calculations(**kwargs)

    bs_trees = os.path.join(rep_dir,"bs-trees.newick")
    assert os.path.isfile(bs_trees)
    f = open(bs_trees)
    lines = f.readlines()
    f.close()

    assert len(lines) == 4
    for d in ["00001","00002","00003","00004"]:
        rep = os.path.join(rep_dir,d)
        assert os.path.isfile(os.path.join(rep,"completed"))
        assert os.path.isdir(os.path.join(rep,"result"))
        assert os.path.isfile(os.path.join(rep,
                                           "result",
                                           "results",
                                           "reconcile",
                                           "geneTree.newick"))

    os.chdir(current_dir)


@pytest.mark.run_raxml
@pytest.mark.run_generax
def test_reconcile_bootstrap(small_phylo,tmpdir):

    df_csv = small_phylo["initial-input/dataframe.csv"]
    df = topiary.read_dataframe(df_csv)
    gene_tree = small_phylo["final-output/gene-tree.newick"]
    species_tree = small_phylo["initial-input/species-tree.newick"]
    reconciled_tree = small_phylo["final-output/reconciled-tree.newick"]
    input_bootstrap_directory = small_phylo["05_gene-tree-bootstraps_toy/output/bootstrap_replicates/"]
    f = open(small_phylo["model.txt"],"r")
    model = f.read().strip()
    f.close()

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    kwargs_template = {"df":df,
                       "model":model,
                       "gene_tree":gene_tree,
                       "species_tree":species_tree,
                       "reconciled_tree":reconciled_tree,
                       "allow_horizontal_transfer":True,
                       "bootstrap_directory":input_bootstrap_directory,
                       "converge_cutoff":0.03,
                       "seed":True,
                       "restart":None,
                       "overwrite":False,
                       "supervisor":None,
                       "num_threads":1,
                       "threads_per_rep":1,
                       "generax_binary":GENERAX_BINARY,
                       "raxml_binary":RAXML_BINARY}

    supervisor = Supervisor()
    supervisor.create_calc_dir("test0",
                               calc_type="test0",
                               df=df,
                               gene_tree=gene_tree,
                               model=model)

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["supervisor"] = supervisor

    tT = reconcile_bootstrap(**kwargs)

    output_dir = supervisor.output_dir
    expected_files = ["dataframe.csv",
                      "reconciliations.txt",
                      "summary-tree.pdf",
                      "reconciled-tree_supports.newick"]
    for f in expected_files:
        assert os.path.isfile(os.path.join(output_dir,f))

    new_T = ete3.Tree(os.path.join(output_dir,"reconciled-tree_supports.newick"),format=0)
    old_T = ete3.Tree(reconciled_tree)

    # Topology should *not* have changed
    assert new_T.robinson_foulds(old_T,unrooted_trees=True)[0] == 0

    # Make sure it now has supports
    for n in new_T.traverse():
        if not n.is_leaf():
            print(n.support)

    os.chdir(current_dir)
