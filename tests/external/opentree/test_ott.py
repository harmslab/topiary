
import pytest
import topiary
from topiary import opentree

import pandas as pd
import numpy as np

import re

def test_get_ott_id(test_dataframes):

    df = test_dataframes["good-df"]
    out_df = opentree.ott.get_ott_id(df)
    assert out_df is not df

    tmp_df = df.drop(columns="ott")
    out_df = opentree.ott.get_ott_id(df)
    assert np.array_equal(out_df.loc[:,"ott"],
                          df.loc[:,"ott"])
