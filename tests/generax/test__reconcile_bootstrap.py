import pytest
import topiary

from topiary.generax._reconcile_bootstrap import _combine_results
from topiary.generax._reconcile_bootstrap import _create_bootstrap_dirs
from topiary.generax._reconcile_bootstrap import reconcile_bootstrap
from topiary.generax._generax import GENERAX_BINARY

import ete3

import os
import glob
import shutil

def test__create_bootstrap_dirs(generax_data,tmpdir):

    input_dir = os.path.join(generax_data["toy-input"],"toy-bootstrap")
    output = os.path.join(tmpdir,"toy-reconcile-bootstrap")

    _create_bootstrap_dirs(previous_dir=input_dir,
                           output=output)

    # Make sure that directories are made with correct files
    expected_files = set(["alignment.phy",
                          "control.txt",
                          "gene_tree.newick",
                          "mapping.link",
                          "run_generax.sh",
                          "species_tree.newick"])
    expected = os.path.join(output,"replicates")
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
        alignments.append(_read_file(os.path.join(input_dir,"output","bootstrap_replicates",f"bsmsa_000{i+1}.phy")))

    trees = []
    with open(os.path.join(input_dir,"output","bootstrap_replicates","bootstraps.newick")) as f:
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

def test__combine_results():

    # hard to test this without saving intermediate outputs and checking. Mostly
    # checked by final reconcilation. Deprioritizing test
    pass


def test_reconcile_bootstrap(generax_data,tmpdir):

    input_dir = os.path.join(generax_data["toy-input"],"toy-bootstrap")
    output = os.path.join(tmpdir,"toy-reconcile-bootstrap-full")
    if os.path.exists(output):
        shutil.rmtree(output)

    kwargs = {"previous_dir":input_dir,
              "df":None,
              "model":None,
              "tree_file":None,
              "species_tree_file":None,
              "allow_horizontal_transfer":True,
              "output":output,
              "overwrite":False,
              "num_threads":2,
              "generax_binary":GENERAX_BINARY}

    tT = reconcile_bootstrap(**kwargs)

    output_dir = os.path.join(output,"output")
    expected_files = ["dataframe.csv","reconciliations.txt",
                      "run_parameters.json","summary-tree.pdf","tree.newick",
                      "tree_supports.newick"]
    for f in expected_files:
        assert os.path.isfile(os.path.join(output_dir,f))

    new_T = ete3.Tree(os.path.join(output_dir,"tree_supports.newick"),format=0)
    old_T = ete3.Tree(os.path.join(input_dir,"output","tree.newick"))

    # Topology should *not* have changed
    assert new_T.robinson_foulds(old_T,unrooted_trees=True)[0] == 0

    # Make sure it now has supports
    for n in new_T.traverse():
        if not n.is_leaf():
            print(n.support)
