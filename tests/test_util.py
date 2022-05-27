
import pytest
import topiary
from topiary import util

import pandas as pd
import numpy as np

import re

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

    # Check alignment
    df = util.check_topiary_dataframe(test_dataframes["good-df_with-good-alignment"])
    df["alignment"]

    # Send in some bad alignments
    bad_align_keys = ["bad-alignment-length1","bad-alignment-length2"]
    for k in bad_ott_keys:
        with pytest.raises(ValueError):
            df = util.check_topiary_dataframe(test_dataframes[k])

    # Check for gap-only column deletion
    input_df = test_dataframes["good-df_with-gap-only-col-alignment"]
    checked_df = util.check_topiary_dataframe(input_df)
    assert input_df["alignment"].iloc[0] == "MLPFLFF---"
    assert checked_df["alignment"].iloc[0] == "MLPFLFF--"
    assert input_df["alignment"].iloc[-1] == "MLPFLFF-TL"
    assert checked_df["alignment"].iloc[-1] == "MLPFLFFTL"

def test__compile_paralog_patterns():

    # Send in bad paralog_patterns data types (should be dict)
    bad_inputs = [1,-1,1.5,False,pd.DataFrame]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            patterns = util._compile_paralog_patterns(paralog_patterns=b)

    # Send bad paralog_patterns keys (should be bad)
    bad_inputs = [(1,2),-1,False]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            patterns = util._compile_paralog_patterns(paralog_patterns=b)

    # Send bad paralog_patterns values
    bad_inputs = [(1,2),-1,False]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            patterns = util._compile_paralog_patterns(paralog_patterns=b)

    # Various paralog_patterns calls hould work
    patterns = util._compile_paralog_patterns(paralog_patterns={"test":"this"})
    assert len(patterns) == 1
    assert patterns[0][0].search("string matches this") is not None
    assert patterns[0][0].search("string does not match") is None
    assert patterns[0][1] == "test"

    patterns = util._compile_paralog_patterns(paralog_patterns={"test":["this","other"]})
    assert len(patterns) == 1
    assert patterns[0][0].search("string matches this") is not None
    assert patterns[0][0].search("string matches other") is not None
    assert patterns[0][0].search("string does not match") is None
    assert patterns[0][1] == "test"

    patterns = util._compile_paralog_patterns(paralog_patterns={"test":("this","other")})
    assert len(patterns) == 1
    assert patterns[0][0].search("string matches this") is not None
    assert patterns[0][0].search("string matches other") is not None
    assert patterns[0][0].search("string does not match") is None
    assert patterns[0][1] == "test"

    patterns = util._compile_paralog_patterns(paralog_patterns={"test":(re.compile("this"),"other")})
    assert len(patterns) == 1
    assert patterns[0][0].search("string matches this") is not None
    assert patterns[0][0].search("string matches other") is not None
    assert patterns[0][0].search("string does not match") is None
    assert patterns[0][1] == "test"

    patterns = util._compile_paralog_patterns(paralog_patterns={"test":"this",
                                                                "thing":"other"})
    assert len(patterns) == 2
    assert patterns[0][0].search("string matches this") is not None
    assert patterns[0][0].search("string matches other") is None
    assert patterns[0][1] == "test"
    assert patterns[1][0].search("string matches other") is not None
    assert patterns[1][0].search("string does not match") is None
    assert patterns[1][1] == "thing"

    patterns = util._compile_paralog_patterns(paralog_patterns={"test":"this"})
    assert patterns[0][0].flags == re.compile("a",flags=re.IGNORECASE).flags

    patterns = util._compile_paralog_patterns(paralog_patterns={"test":"this"},
                                              ignorecase=True)
    assert patterns[0][0].flags == re.compile("a",flags=re.IGNORECASE).flags

    patterns = util._compile_paralog_patterns(paralog_patterns={"test":"this"},
                                              ignorecase=False)
    assert patterns[0][0].flags == re.compile("a").flags


def test_create_nicknames(test_dataframes):

    df = test_dataframes["good-df"]

    # Make sure it runs on a copy
    out_df = util.create_nicknames(df,paralog_patterns={})
    assert out_df is not df

    # Check source column check
    with pytest.raises(ValueError):
        out_df = util.create_nicknames(df,paralog_patterns={},source_column="not_a_column")

    # Check trying to overwrite reserved column
    with pytest.raises(ValueError):
        out_df = util.create_nicknames(df,paralog_patterns={},output_column="ott")

    # Check trying to overwrite reserved column -- should still throw error
    # even if we try to force overwrite
    with pytest.raises(ValueError):
        out_df = util.create_nicknames(df,paralog_patterns={},output_column="ott",overwrite_output=True)

    # Throw error to avoid overwrite
    with pytest.raises(ValueError):
        out_df = util.create_nicknames(df,paralog_patterns={},output_column="isoform")

    # Force overwrite
    out_df = util.create_nicknames(df,paralog_patterns={},output_column="isoform",overwrite_output=True)

    # Send in bad paralog_patterns separator (should be string)
    bad_inputs = [1,-1,1.5,False,pd.DataFrame]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            out_df = util.create_nicknames(df,paralog_patterns={},separator=b)

    # Should work without error
    out_df = util.create_nicknames(df,paralog_patterns={},separator="a string")

    # Send in bad unassigned  (should be string)
    bad_inputs = [1,-1,1.5,False,pd.DataFrame]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            out_df = util.create_nicknames(df,paralog_patterns={},unassigned_name=b)

    # Should work without error
    out_df = util.create_nicknames(df,paralog_patterns={},unassigned_name="a string")

    # Send in bad paralog_patterns data types (should be dict)
    bad_inputs = [1,-1,1.5,False,pd.DataFrame]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            out_df = util.create_nicknames(df,paralog_patterns=b)

    # Send bad paralog_patterns keys (should be bad)
    bad_inputs = [(1,2),-1,False]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            out_df = util.create_nicknames(df,paralog_patterns={b:"test"})

    # Send bad paralog_patterns values
    bad_inputs = [(1,2),-1,False]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            out_df = util.create_nicknames(df,paralog_patterns={"test":["this",b]})

    # Various paralog_patterns calls hould work
    out_df = util.create_nicknames(df,paralog_patterns={"test":"this"})
    out_df = util.create_nicknames(df,paralog_patterns={"test":["this","this"]})
    out_df = util.create_nicknames(df,paralog_patterns={"test":("this","this")})
    out_df = util.create_nicknames(df,paralog_patterns={"test":(re.compile("this"),"this")})

    # Check to make nicknaming is doing what we expect
    test_df = df.copy()
    test_df.loc[:,"name"] = ["rocking","out","in","the","usa"]
    paralog_patterns = {"fixed":("rock","out"),
                        "junk":("the",re.compile("usa"))}
    out_df = util.create_nicknames(test_df,output_column="test1",paralog_patterns=paralog_patterns)

    assert np.array_equal(np.array(out_df.loc[:,"test1"]),
                          np.array(["fixed","fixed","unassigned","junk","junk"]))

    # Do again, but do some stuff that requires escape characters and more
    # interesting regular expressions
    test_df = df.copy()
    test_df.loc[:,"name"] = ["|rocking","rock","\in","the","usa"]
    paralog_patterns = {"fixed":("|rock","out"),
                        "junk":("the",re.compile("us."))}
    out_df = util.create_nicknames(test_df,output_column="test1",paralog_patterns=paralog_patterns)

    assert np.array_equal(np.array(out_df.loc[:,"test1"]),
                          np.array(["fixed","unassigned","unassigned","junk","junk"]))

    # Make sure ignorecase is done correctly (this should work because ignorecase
    # defaults to True)
    test_df = df.copy()
    test_df.loc[:,"name"] = ["rocking","OUT","in","the","usa"]
    paralog_patterns = {"fixed":("rock","out"),
                        "junk":("the","usa")}
    out_df = util.create_nicknames(test_df,output_column="test1",paralog_patterns=paralog_patterns)

    assert np.array_equal(np.array(out_df.loc[:,"test1"]),
                          np.array(["fixed","fixed","unassigned","junk","junk"]))

    # Make sure ignorecase is done correctly (this should not replace "OUT"
    # because ignorecase is False
    test_df = df.copy()
    test_df.loc[:,"name"] = ["rocking","OUT","in","the","usa"]
    paralog_patterns = {"fixed":("rock","out"),
                        "junk":("the","usa")}
    out_df = util.create_nicknames(test_df,output_column="test1",paralog_patterns=paralog_patterns,
                                   ignorecase=False)

    assert np.array_equal(np.array(out_df.loc[:,"test1"]),
                          np.array(["fixed","unassigned","unassigned","junk","junk"]))

    # Make sure we can control the source column
    test_df = df.copy()
    paralog_patterns = {"froggy":"Hylobates"}
    out_df = util.create_nicknames(test_df,source_column="species",paralog_patterns=paralog_patterns)
    assert out_df["nickname"].iloc[0] == "froggy"
