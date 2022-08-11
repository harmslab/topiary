
import pytest
import topiary

from topiary.generax._reconcile_no_bootstrap import reconcile_no_bootstrap
from topiary.generax._generax import GENERAX_BINARY

import ete3

import os
import shutil

def test_reconcile_no_bootstrap(generax_data,tmpdir):

    input_dir = os.path.join(generax_data["toy-input"],"toy-ml")
    output = os.path.join(tmpdir,"toy-reconcile-full")
    if os.path.exists(output):
        shutil.rmtree(output)

    # Send in a wrong tree -- should be fixed by correct generax run!
    kwargs = {"previous_dir":input_dir,
              "df":None,
              "model":None,
              "tree_file":os.path.join(input_dir,"output","tree_wrong.newick"),
              "species_tree_file":None,
              "allow_horizontal_transfer":True,
              "output":output,
              "overwrite":False,
              "num_threads":2,
              "generax_binary":GENERAX_BINARY}

    tT = reconcile_no_bootstrap(**kwargs)

    output_dir = os.path.join(output,"output")
    expected_files = ["dataframe.csv",
                      "run_parameters.json","summary-tree.pdf","tree.newick"]
    for f in expected_files:
        assert os.path.isfile(os.path.join(output_dir,f))

    assert os.path.isdir(os.path.join(output_dir,"reconcilations"))

    new_T = ete3.Tree(os.path.join(output_dir,"tree.newick"),format=0)
    old_T_wrong = ete3.Tree(os.path.join(input_dir,"output","tree_wrong.newick"))
    old_T = ete3.Tree(os.path.join(output_dir,"tree.newick"),format=0)

    # Toplogy should have changed to be correct
    assert new_T.robinson_foulds(old_T_wrong,unrooted_trees=True)[0] != 0

    # Correct toplogy
    assert new_T.robinson_foulds(old_T,unrooted_trees=True)[0] == 0
