import pytest
import topiary

from topiary.generax._reconcile_bootstrap import _create_bootstrap_dirs
from topiary.generax._reconcile_bootstrap import _run_bootstrap_calculations
from topiary.generax._reconcile_bootstrap import _combine_bootstrap_calculations
from topiary.generax._reconcile_bootstrap import reconcile_bootstrap
from topiary.generax._generax import GENERAX_BINARY
from topiary._private import Supervisor

import ete3

import os
import glob
import shutil
import copy

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test__check_calc_completeness():
    # Status bar runs on it's own thread. Low priority (and pain) to test.
    pass


@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test__create_bootstrap_dirs(generax_data,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    input_dir = os.path.abspath(os.path.join(generax_data["toy-input"],"toy-bootstrap","output"))
    df = topiary.read_dataframe(os.path.join(input_dir,"dataframe.csv"))
    model = "JTT"
    tree_file = os.path.join(input_dir,"tree.newick")
    species_tree_file = os.path.join(input_dir,"species_tree.newick")
    bootstrap_directory = os.path.join(input_dir,"bootstrap_replicates")

    kwargs_template = {"df":df,
                       "model":model,
                       "tree_file":tree_file,
                       "species_tree_file":species_tree_file,
                       "allow_horizontal_transfer":True,
                       "bootstrap_directory":bootstrap_directory,
                       "overwrite":False,
                       "generax_binary":GENERAX_BINARY}

    os.mkdir("test0")
    os.chdir("test0")
    kwargs = copy.deepcopy(kwargs_template)

    _create_bootstrap_dirs(**kwargs)

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
    with open(os.path.join(input_dir,"bootstrap_replicates","bootstraps.newick")) as f:
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

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test__run_bootstrap_calculations(generax_data,tmpdir):

    # This basically makes sure that a calculation is run in every directory.
    # It does not check quality/correctness of output

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    input_dir = os.path.abspath(os.path.join(generax_data["toy-input"],"toy-bootstrap"))
    replicate_dir = os.path.join(input_dir,"generax-bs-replicates_pre-run")

    shutil.copytree(replicate_dir,"replicates")

    _run_bootstrap_calculations(replicate_dir="replicates",
                                num_threads=1)

    for d in ["00001","00002","00003","00004"]:
        rep = os.path.join("replicates",d)
        ran = glob.glob(os.path.join(rep,"completed_run-*"))
        assert len(ran) == 1
        assert os.path.isdir(os.path.join(rep,"result"))
        assert os.path.isfile(os.path.join(rep,
                                           "result",
                                           "results",
                                           "reconcile",
                                           "geneTree.newick"))

    os.chdir(current_dir)

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test__combine_bootstrap_calculations(generax_data,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    input_dir = os.path.abspath(os.path.join(generax_data["toy-input"],"toy-bootstrap"))
    replicate_dir = os.path.join(input_dir,"generax-bs-replicates_post-run")
    tree_file = os.path.join(input_dir,"output","tree.newick")

    shutil.copytree(replicate_dir,"replicates")

    converged = _combine_bootstrap_calculations("replicates",tree_file)

    # Make sure we made a bs-trees.newick file with 4 different lines
    f = open("bs-trees.newick")
    lines = f.readlines()
    f.close()
    assert len(lines) == 4
    for i in range(3):
        assert lines[i] != lines[3]

    # Make sure that we loaded supports on that are different from each other
    supports = []
    T = ete3.Tree("tree_supports.newick")
    for n in T.traverse():
        if not n.is_leaf():
            supports.append(n.support)
    assert len(set(supports)) > 1

    assert isinstance(converged,bool)

    os.chdir(current_dir)


@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_reconcile_bootstrap(generax_data,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    input_dir = os.path.abspath(os.path.join(generax_data["toy-input"],"toy-bootstrap","output"))
    df = topiary.read_dataframe(os.path.join(input_dir,"dataframe.csv"))
    model = "JTT"
    tree_file = os.path.join(input_dir,"tree.newick")
    species_tree_file = os.path.join(input_dir,"species_tree.newick")
    bootstrap_directory = os.path.join(input_dir,"bootstrap_replicates")

    kwargs_template = {"df":df,
                       "model":model,
                       "tree_file":tree_file,
                       "species_tree_file":species_tree_file,
                       "allow_horizontal_transfer":True,
                       "bootstrap_directory":bootstrap_directory,
                       "overwrite":False,
                       "supervisor":None,
                       "num_threads":1,
                       "generax_binary":GENERAX_BINARY}

    supervisor = Supervisor()
    supervisor.create_calc_dir("test0",
                               calc_type="test0",
                               df=df,
                               tree=tree_file,
                               model=model)

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["supervisor"] = supervisor

    tT = reconcile_bootstrap(**kwargs)

    output_dir = supervisor.output_dir
    expected_files = ["dataframe.csv",
                      "reconciliations.txt",
                      "summary-tree.pdf",
                      "tree_supports.newick"]
    for f in expected_files:
        assert os.path.isfile(os.path.join(output_dir,f))

    new_T = ete3.Tree(os.path.join(output_dir,"tree_supports.newick"),format=0)
    old_T = ete3.Tree(os.path.join(input_dir,"tree.newick"))

    # Topology should *not* have changed
    assert new_T.robinson_foulds(old_T,unrooted_trees=True)[0] == 0

    # Make sure it now has supports
    for n in new_T.traverse():
        if not n.is_leaf():
            print(n.support)

    os.chdir(current_dir)
