
import pytest

from topiary.raxml.convergence import _parse_convergence_file
from topiary.raxml.convergence import check_convergence

import numpy as np
import pandas as pd

import os
import shutil

def test__parse_convergence_file(raxml_output):

    c, df = _parse_convergence_file(raxml_output["bs-trees.log_converged"])
    assert c
    assert len(df) == 6
    assert np.array_equal(df.columns,["trees","avg_wrf","avg_wrf_pct",
                                      "perms_below_cutoff","converged"])

    c, df = _parse_convergence_file(raxml_output["bs-trees.log_not-converged"])
    assert not c
    assert len(df) == 20
    assert np.array_equal(df.columns,["trees","avg_wrf","avg_wrf_pct",
                                      "perms_below_cutoff","converged"])

@pytest.mark.run_raxml
def test_check_convergence(tiny_phylo,tmpdir):

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    shutil.copy(tiny_phylo["05_reconcile-bootstraps_toy/working/bs-trees.newick"],
                "bs-trees.newick")

    c, df = check_convergence("bs-trees.newick")
    assert not c
    assert issubclass(type(df),pd.DataFrame)
    assert not os.path.exists("tmp_test_conv")

    # Make sure we can write out directory
    c, df = check_convergence("bs-trees.newick",calc_dir="save-tmp")
    assert not c
    assert issubclass(type(df),pd.DataFrame)
    assert os.path.exists("save-tmp")

    with open("save-tmp/bs-trees.newick.raxml.log") as f:
        for line in f:
            if line.startswith("# trees"):
                cutoff = line.split()[12].strip()
                assert cutoff == "3.00"
                break
    shutil.rmtree("save-tmp")

    # Make sure the cutoff is being set
    c, df = check_convergence("bs-trees.newick",
                              calc_dir="save-tmp",
                              converge_cutoff=0.05)
    assert not c
    assert issubclass(type(df),pd.DataFrame)
    assert os.path.exists("save-tmp")

    with open("save-tmp/bs-trees.newick.raxml.log") as f:
        for line in f:
            if line.startswith("# trees"):
                cutoff = line.split()[12].strip()
                assert cutoff == "5.00"
                break
    shutil.rmtree("save-tmp")

    # Make sure the seed is being set
    c, df = check_convergence("bs-trees.newick",
                              calc_dir="save-tmp",
                              seed=12345)
    assert not c
    assert issubclass(type(df),pd.DataFrame)
    assert os.path.exists("save-tmp")

    with open("save-tmp/bs-trees.newick.raxml.log") as f:
        for line in f:
            if line.startswith("random"):
                seed = line.split()[2].strip()
                assert seed == "12345"
                break
    shutil.rmtree("save-tmp")


    os.chdir(current_dir)
