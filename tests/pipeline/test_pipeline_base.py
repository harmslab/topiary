
import pytest

import topiary

import numpy as np
import pandas as pd

import os

def _a_dummy_function(df,stupid,does):
    """
    Dummy function for testing _run_and_print
    """

    df["species"] = stupid
    df["accession"] = does

    return df

def test__run_and_print(test_dataframes,tmpdir):
    """
    Test internal _run_and_print. Because internal just check simple functionality,
    not arg checking.
    """

    from topiary.pipeline.base import _run_and_print as rap

    input_df = test_dataframes["good-df"].copy()

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    df, i = rap(function=_a_dummy_function,
                kwargs={"df":input_df,"stupid":"Bob Dole","does":5},
                step_counter=0,
                out_file_string="ya-basic",
                human_string="message")

    out_file = "00_ya-basic-dataframe.csv"
    file_df = pd.read_csv(out_file)
    os.remove(out_file)

    expected_species = ["Bob Dole" for _ in range(len(input_df))]
    expected_acc = [5 for _ in range(len(input_df))]

    assert np.array_equal(df["species"],expected_species)
    assert np.array_equal(file_df["species"],expected_species)
    assert np.array_equal(df["accession"],expected_acc)
    assert np.array_equal(file_df["accession"],expected_acc)
    assert i == 1

    os.chdir(current_dir)
