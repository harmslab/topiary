
import pytest
import topiary
from topiary import util

import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------- #
# Test __init__
# ---------------------------------------------------------------------------- #

def test_create_pipeline_dict():

    pipe_dict = topiary.util.create_pipeline_dict()

    key_list = ["name","species","sequence","uid","keep"]

    assert set(pipe_dict.keys()) == set(key_list)

def test_grab_line_meta_data(ncbi_lines,ncbi_lines_parsed):

    for i, line in enumerate(ncbi_lines):
        line_dict = topiary.util.grab_line_meta_data(line)

        out = []
        for k in ncbi_lines_parsed[i]:
            try:
                out.append(line_dict[k] == ncbi_lines_parsed[i][k])
            except KeyError:
                pass

        assert sum(out) == len(ncbi_lines_parsed[i])

#def pretty_to_uid(df,to_convert,out_file=None,overwrite=False):
#def uid_to_pretty(df,to_convert,out_file=None,overwrite=False):
#def load_tree(tree,fmt=None):

def test_check_topiary_dataframe(test_dataframes):
    """
    Test check for topiary dataframe.
    """

    # Make sure reads good dataframe without mangling
    good_df = test_dataframes["good-df"]
    df = util.check_topiary_dataframe(test_dataframes["good-df"])

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
            util.check_topiary_dataframe(b)

    # Make sure it drops empty lines
    df = util.check_topiary_dataframe(test_dataframes["good-test_blank-lines"])
    assert len(df) == 5

    # Make sure it properly looks for required columns
    required_col = ["species","name","sequence"]
    for c in required_col:
        bad_df = good_df.copy()
        bad_df = bad_df.drop(columns=[c])
        with pytest.raises(ValueError):
            util.check_topiary_dataframe(bad_df)

    # Check keep
    df = util.check_topiary_dataframe(test_dataframes["good-test-keep-parse"])
    assert df.keep.dtypes is np.dtype(bool)
    assert np.array_equal(df.keep,np.array([True,False,
                                            True,False,
                                            True,False]))

    df = util.check_topiary_dataframe(test_dataframes["good-test-keep-parse_number"])
    assert df.keep.dtypes is np.dtype(bool)
    assert np.array_equal(df.keep,np.array([True,False,
                                            True,False]))


    # Check add uid column.
    df = util.check_topiary_dataframe(test_dataframes["no-uid"])
    assert len(np.unique(df.uid)) == len(df)

    # Check replace bad uid.
    df = util.check_topiary_dataframe(test_dataframes["bad-uid"])
    for d in df.uid:
        assert type(d) is str
        assert len(d) == 10

    # Check make uid unique.
    df = util.check_topiary_dataframe(test_dataframes["duplicate-uid"])
    assert len(np.unique(df.uid)) == len(df)

    # Check ott
    bad_ott_keys = ["bad-ott1","bad-ott2","bad-ott3"]
    for k in bad_ott_keys:
        with pytest.raises(ValueError):
            df = util.check_topiary_dataframe(test_dataframes[k])

def test_create_nicknames(test_dataframes):

    df = test_dataframes["good-df"]

    # Make sure it runs on a copy
    out_df = util.create_nicknames(df)
    assert out_df is not df

    with pytest.raises(ValueError):
        out_df = util.create_nicknames(df,source_column="not_a_column")

    with pytest.raises(ValueError):
        out_df = util.create_nicknames(df,output_column="ott")

    # Throw error to avoid overwrite
    with pytest.raises(ValueError):
        out_df = util.create_nicknames(df,output_column="isoform")

    # Force overwrite
    out_df = util.create_nicknames(df,output_column="isoform",overwrite_output=True)

    # Send in bad alias data types
    bad_inputs = [1,-1,1.5,False,pd.DataFrame]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            out_df = util.create_nicknames(df,aliases=b)

    # Send bad aliases keys
    bad_inputs = [(1,2),-1,False]
    with pytest.raises(ValueError):
        out_df = util.create_nicknames(df,aliases={b:"test"})

    # Send bad aliases values
    bad_inputs = [(1,2),-1,False]
    with pytest.raises(ValueError):
        out_df = util.create_nicknames(df,aliases={"test":["this",b]})

    # Various alias calls hould work
    out_df = util.create_nicknames(df,aliases={"test":"this"})
    out_df = util.create_nicknames(df,aliases={"test":["this","this"]})
    out_df = util.create_nicknames(df,aliases={"test":("this","this")})

    # Check to make sure replacement happening
    test_df = df.copy()
    test_df.loc[:,"name"] = ["rocking","out","in","the","usa"]
    aliases = {"fixed":("rock","out"),
               "junk":("the","usa")}
    out_df = util.create_nicknames(test_df,output_column="test1",aliases=aliases)

    assert np.array_equal(np.array(out_df.loc[:,"test1"]),
                          np.array(["fixeding","fixed","in","junk","junk"]))

    # Check auto_remove_common functionality
    test_df = df.copy()
    test_df.loc[:,"name"] = "protein powder"
    out_df = util.create_nicknames(test_df,output_column="test1")
    assert out_df["test1"].iloc[0] == "powder"

    out_df = util.create_nicknames(test_df,output_column="test1",auto_remove_common=False)
    assert out_df["test1"].iloc[0] == "protein powder"

    test_df = df.copy()
    aliases = {"froggy":"Hylobates"}
    out_df = util.create_nicknames(test_df,source_column="species",aliases=aliases)
    assert out_df["nickname"].iloc[0] == "froggy moloch"

def test_get_ott_id(test_dataframes):

    df = test_dataframes["good-df"]
    out_df = util.get_ott_id(df)
    assert out_df is not df

    tmp_df = df.drop(columns="ott")
    out_df = util.get_ott_id(df)
    assert np.array_equal(out_df.loc[:,"ott"],
                          df.loc[:,"ott"])
