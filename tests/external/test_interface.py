import pytest

import topiary
import topiary.external.interface as interface
import numpy as np
import pandas as pd

import warnings, os

def test_launch(tmpdir):
    """
    Test the launcher.
    """

    pass

def test_read_previous_run_dir(run_directories):

    bad_inputs = [1,None,{"test":1},(),[]]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            interface.read_previous_run_dir(b)

    with pytest.raises(ValueError):
        interface.read_previous_run_dir("not_really_a_directory")

    # Pass something too deep into good run directory
    with pytest.raises(ValueError):
        interface.read_previous_run_dir(os.path.join(run_directories["ml-tree"],
                                              "output"))

    # no json; should die
    with pytest.raises(ValueError):
        interface.read_previous_run_dir(run_directories["ml-tree_no-json"])

    # bad dataframe; should die (make sure it's checking df on way in)
    with pytest.raises(ValueError):
        interface.read_previous_run_dir(run_directories["ml-tree_bad-df"])

    # no dataframe; should die
    with pytest.raises(ValueError):
        interface.read_previous_run_dir(run_directories["ml-tree_no-df"])

    # no tree; should work fine
    interface.read_previous_run_dir(run_directories["ml-tree_no-tree"])

    # Make sure the function is properly reading run directories
    for t in [("ml-tree","ml_tree"),
              ("ancestors","ancestors"),
              ("find-model","find_best_model"),
              ("reconciled","reconciliation")]:

        key = t[0]
        calc_type = t[1]

        out = interface.read_previous_run_dir(run_directories[key])
        assert out["model"] == "JTT"
        assert out["calc_type"] == calc_type
        assert type(out["cmd"]) in [type("string"),type(None)]
        try:
            tf = out["tree_file"]
            assert os.path.exists(tf)
        except KeyError:
            pass
        # This will throw an error if the dataframe read is bad
        topiary.util.check_topiary_dataframe(out["df"])
