import pytest

import topiary
from topiary.check import check_topiary_dataframe

import numpy as np
import pandas as pd

def test_check_topiary_dataframe(test_dataframes):
    """
    Test check for topiary dataframe.
    """

    # Make sure reads good dataframe without mangling
    good_df = test_dataframes["good-df"]
    df = check_topiary_dataframe(test_dataframes["good-df"])

    # Sort by column name
    orig_columns = list(good_df.columns)
    orig_columns.sort()
    new_columns = list(df.columns)
    new_columns.sort()

    assert np.sum(np.asarray(good_df.loc[:,orig_columns] == df.loc[:,new_columns]) == False) == 0

    # make sure it properly deals with all sorts of wacky input
    bad_inputs = [1,-1,1.5,None,False,pd.DataFrame]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            check_topiary_dataframe(b)

    # Make sure it drops empty lines
    df = check_topiary_dataframe(test_dataframes["good-test_blank-lines"])
    assert len(df) == 5

    # Make sure it properly looks for required columns
    required_col = ["species","name","sequence"]
    for c in required_col:
        bad_df = good_df.copy()
        bad_df = bad_df.drop(columns=[c])
        with pytest.raises(ValueError):
            check_topiary_dataframe(bad_df)

    # Make sure it splices strings with \n and leading/trailing " "
    test_df = good_df.copy()
    test_df.loc[:,"sequence"] = ["TESTTHIS","TEST THIS","TEST\nTHIS"," TESTTHIS "," TEST\n THIS "]
    test_df = test_df.drop(columns=["length"])
    df = check_topiary_dataframe(test_df)
    assert np.array_equal(df.sequence,["TESTTHIS" for _ in range(5)])

    # Check keep
    df = check_topiary_dataframe(test_dataframes["good-test-keep-parse"])
    assert df.keep.dtypes is np.dtype(bool)
    assert np.array_equal(df.keep,np.array([True,False,
                                            True,False,
                                            True,False]))

    df = check_topiary_dataframe(test_dataframes["good-test-keep-parse_number"])
    assert df.keep.dtypes is np.dtype(bool)
    assert np.array_equal(df.keep,np.array([True,False,
                                            True,False]))


    # Check add uid column.
    df = check_topiary_dataframe(test_dataframes["no-uid"])
    assert len(np.unique(df.uid)) == len(df)

    # Check replace bad uid.
    df = check_topiary_dataframe(test_dataframes["bad-uid"])
    for d in df.uid:
        assert type(d) is str
        assert len(d) == 10

    # Check make uid unique.
    df = check_topiary_dataframe(test_dataframes["duplicate-uid"])
    assert len(np.unique(df.uid)) == len(df)

    # Check ott
    bad_ott_keys = ["bad-ott1","bad-ott2","bad-ott3"]
    for k in bad_ott_keys:
        with pytest.raises(ValueError):
            df = check_topiary_dataframe(test_dataframes[k])

    # Check alignment
    df = check_topiary_dataframe(test_dataframes["good-df_with-good-alignment"])
    df["alignment"]

    # Send in some bad alignments
    bad_align_keys = ["bad-alignment-length1","bad-alignment-length2"]
    for k in bad_ott_keys:
        with pytest.raises(ValueError):
            df = check_topiary_dataframe(test_dataframes[k])

    # Check for gap-only column deletion
    input_df = test_dataframes["good-df_with-gap-only-col-alignment"]
    checked_df = check_topiary_dataframe(input_df)
    assert input_df["alignment"].iloc[0] == "MLPFLFF---"
    assert checked_df["alignment"].iloc[0] == "MLPFLFF--"
    assert input_df["alignment"].iloc[-1] == "MLPFLFF-TL"
    assert checked_df["alignment"].iloc[-1] == "MLPFLFFTL"
