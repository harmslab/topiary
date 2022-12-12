
import pytest
import topiary
from topiary.opentree.ott import get_df_ott

import pandas as pd
import numpy as np

import re

def test_get_df_ott(test_dataframes):


    df = test_dataframes["good-df"]
    out_df = get_df_ott(df)
    assert out_df is not df

    tmp_df = df.drop(columns="ott")
    out_df = get_df_ott(tmp_df)
    assert np.array_equal(out_df.loc[:,"ott"],df.loc[:,"ott"])
    assert np.array_equal(out_df.loc[:,"species"],df.loc[:,"orig_species"])
    assert np.array_equal(np.ones(len(tmp_df),dtype=bool),out_df.keep)

    # make sure that it handles bad species names gracefully, setting keep to
    # False
    tmp_df = df.drop(columns="ott")
    tmp_df.loc[tmp_df.index[0],"species"] = "Not a species"
    out_df = get_df_ott(tmp_df)
    assert pd.isnull(out_df.loc[out_df.index[0],"ott"])
    assert out_df.loc[out_df.index[0],"species"] == "Not a species"
    expected_keep = np.ones(len(out_df),dtype=bool)
    expected_keep[0] = False
    assert np.array_equal(expected_keep,out_df.keep)

    # make sure that it handles all bad species names gracefully, setting keep to
    # False
    tmp_df = df.drop(columns="ott")
    tmp_df.loc[:,"species"] = "Not a species"
    out_df = get_df_ott(tmp_df)
    assert np.sum(pd.isnull(out_df.loc[:,"ott"])) == len(tmp_df)
    assert np.sum(out_df.loc[:,"species"] == "Not a species") == len(tmp_df)
    expected_keep = np.zeros(len(out_df),dtype=bool)
    assert np.array_equal(expected_keep,out_df.keep)

    tmp_df = df.drop(columns=["species","ott"])
    with pytest.raises(ValueError):
        out_df = get_df_ott(tmp_df)

    # make sure that it handles all bad species names gracefully, but keeps the
    # bad ott if we request it. 
    tmp_df = df.drop(columns="ott")
    tmp_df.loc[:,"species"] = "Not a species"
    out_df = get_df_ott(tmp_df,keep_anyway=True)
    assert np.sum(pd.isnull(out_df.loc[:,"ott"])) == len(tmp_df)
    assert np.sum(out_df.loc[:,"species"] == "Not a species") == len(tmp_df)
    expected_keep = np.ones(len(out_df),dtype=bool)
    assert np.array_equal(expected_keep,out_df.keep)


