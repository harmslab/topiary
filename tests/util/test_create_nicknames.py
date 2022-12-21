
import pytest
import topiary
from topiary import util

import pandas as pd
import numpy as np

import re


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
                        "junk":(re.compile("the|us."))}
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
