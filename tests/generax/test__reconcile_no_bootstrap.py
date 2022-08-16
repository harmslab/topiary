
import pytest
import topiary

from topiary.generax._reconcile_no_bootstrap import reconcile_no_bootstrap
from topiary.generax._generax import GENERAX_BINARY
from topiary._private import Supervisor

import ete3

import os
import shutil

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_reconcile_no_bootstrap(generax_data,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    input_dir = os.path.join(generax_data["toy-input"],"toy-ml","output")

    supervisor = Supervisor()
    supervisor.create_calc_dir("test0",
                               calc_type="test0",
                               df=os.path.join(input_dir,"dataframe.csv"),
                               gene_tree=os.path.join(input_dir,"tree_wrong.newick"),
                               species_tree=os.path.join(input_dir,"species_tree.newick"))

    # Send in a wrong tree -- should be fixed by correct generax run
    kwargs = {"supervisor":supervisor,
              "df":supervisor.df,
              "model":"JTT",
              "gene_tree":supervisor.gene_tree,
              "species_tree":os.path.join(input_dir,"species_tree.newick"),
              "allow_horizontal_transfer":True,
              "seed":True,
              "overwrite":False,
              "num_threads":2,
              "generax_binary":GENERAX_BINARY}

    tT = reconcile_no_bootstrap(**kwargs)

    expected_files = ["summary-tree.pdf",
                      "gene-tree.newick",
                      "species-tree.newick",
                      "reconciled-tree.newick",
                      "reconciled-tree_events.newick",
                      "dataframe.csv"]
    for f in expected_files:
        assert os.path.isfile(os.path.join(supervisor.output_dir,f))

    assert os.path.isdir(os.path.join(supervisor.output_dir,"reconcilations"))

    output_T = ete3.Tree(os.path.join(supervisor.output_dir,"reconciled-tree.newick"),format=0)
    input_T = ete3.Tree(os.path.join(supervisor.input_dir,"gene-tree.newick"))
    correct_T = ete3.Tree(os.path.join(input_dir,"tree.newick"),format=0)

    # Toplogy should have changed to be correct
    assert output_T.robinson_foulds(input_T,unrooted_trees=True)[0] != 0

    # Correct toplogy
    assert output_T.robinson_foulds(correct_T,unrooted_trees=True)[0] == 0

    os.chdir(current_dir)
