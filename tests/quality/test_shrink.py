
import pytest

import topiary
from topiary.quality.shrink import shrink_redundant
from topiary.quality.shrink import shrink_aligners
from topiary.quality.shrink import shrink_dataset


import numpy as np
import pandas as pd

import os

def test_shrink_in_species():
    pass

def test_shrink_redundant(for_real_inference):

    df = topiary.read_dataframe(for_real_inference["small-pre-redundancy.csv"])

    # Test with default parameters
    df2 = shrink_redundant(df)
    expected_keep = np.array(df.keep)
    expected_keep[1] = False
    assert np.array_equal(expected_keep,df2.keep)

    # check arg checking
    with pytest.raises(ValueError):
        shrink_redundant(df,paralog_column="not_a_column")

    with pytest.raises(ValueError):
        shrink_redundant(df,weighted_paralog_split="stupid")

    with pytest.raises(ValueError):
        shrink_redundant(df,merge_block_size=-1)

    with pytest.raises(ValueError):
        shrink_redundant(df,redundancy_cutoff=-1)

    # Construct test dataset with only two species
    test_df = df.loc[df.species.isin(["Xenopus laevis","Geotrypetes seraphini"]),:].copy()
    test_df.loc[:,"recip_paralog"] = ["LY96","LY96","LY96","LY86","LY86"]
    test_df.loc[:,"always_keep"] = np.zeros(len(test_df),dtype=bool)
    test_df.loc[:,"key_species"] = np.zeros(len(test_df),dtype=bool)

    # Should whack from LY86 and LY96
    out_df = shrink_redundant(test_df,redundancy_cutoff=0.5)
    assert np.array_equal(out_df.keep,np.array([False,True,True,True,False]))

    # Should whack from LY86 and LY96
    out_df = shrink_redundant(test_df,redundancy_cutoff=0.05)
    assert np.array_equal(out_df.keep,np.array([False,False,True,True,False]))

    # Should keep everything
    test_df_2 = test_df.copy()
    test_df_2.loc[:,"always_keep"] = np.ones(len(test_df),dtype=bool)
    out_df = shrink_redundant(test_df_2,redundancy_cutoff=0.05)
    assert np.array_equal(out_df.keep,np.array([True,True,True,True,True]))

    # Should keep only one now that everything is same paralog
    test_df_2 = test_df.copy()
    test_df_2.loc[:,"recip_paralog"] = ["LY96","LY96","LY96","LY96","LY96"]
    out_df = shrink_redundant(test_df_2,redundancy_cutoff=0.05)
    assert np.array_equal(out_df.keep,np.array([False,False,False,True,False]))

    # Should keep everything because all are different paralogs
    test_df_2 = test_df.copy()
    test_df_2.loc[:,"recip_paralog"] = ["A","B","C","D","E"]
    out_df = shrink_redundant(test_df_2,redundancy_cutoff=0.05)
    assert np.array_equal(out_df.keep,np.array([True,True,True,True,True]))

    # Make a sequence get not kept because it is very short
    test_df_2 = test_df.copy()
    test_df_2.loc[[9],"sequence"] = "MEVW"
    out_df = shrink_redundant(test_df_2,redundancy_cutoff=0.5)
    assert np.array_equal(out_df.keep,np.array([False,False,True,True,False]))

    # All have identical sequence and paralog. Take first because it is the
    # best -- not structure
    test_df_2 = test_df.copy()
    test_df_2.loc[:,"recip_paralog"] = "LY96"
    test_df_2.loc[:,"sequence"] = test_df_2.loc[8,"sequence"]
    test_df_2.loc[:,"structure"] = True
    test_df_2.loc[8,"structure"] = False
    out_df = shrink_redundant(test_df_2,redundancy_cutoff=0.05)
    assert np.array_equal(out_df.keep,np.array([True,False,False,False,False]))

    # All have identical sequence and paralog. Take second because it is not
    # low quality -- rest are.
    test_df_2 = test_df.copy()
    test_df_2.loc[:,"recip_paralog"] = "LY96"
    test_df_2.loc[:,"sequence"] = test_df_2.loc[8,"sequence"]
    test_df_2.loc[:,"structure"] = False
    test_df_2.loc[:,"low_quality"] = True
    test_df_2.loc[9,"low_quality"] = False
    out_df = shrink_redundant(test_df_2,redundancy_cutoff=0.05)
    assert np.array_equal(out_df.keep,np.array([False,True,False,False,False]))

def test_shrink_aligners():
    pass

def test_shrink_dataset(for_real_inference,tmpdir):

    # Skip test of windows because muscle cannot be installed by conda
    if os.name == "nt":
        return

    df = topiary.read_dataframe(for_real_inference["small-pre-redundancy.csv"])
    df.loc[df["recip_paralog"] == "unassigned","recip_paralog"] = "LY96"

    # Runs with default
    out = shrink_dataset(df)

    # no easy way to check that this *works* but at least this should not die
    # when other strings will
    out = shrink_dataset(df,paralog_column="nickname")

    # bad paralog_column values
    bad_paralog_column = ["not_a_column",None,14.2,str]
    for b in bad_paralog_column:
        print(b)
        with pytest.raises(ValueError):
            shrink_dataset(df,paralog_column=b)

    # This should not shrink at all. All sequences should make it fine, as
    # stringency is very low and we allow up to 0.12 seqs/per column. 160
    # columns in alignment, 19 sequences. 160*0.12 --> 19.2
    out = shrink_dataset(df,seqs_per_column=0.12,redundancy_cutoff=1.0,sparse_column_cutoff=1.0)
    assert np.sum(out.keep) == 19

    # Get only six sequence out, as we're only allowing 1/160 seqs per column.
    # There will be the four always_keep
    out = shrink_dataset(df,seqs_per_column=1/160)
    assert np.sum(out.keep) == 4

    tmp_df = df.copy()
    tmp_df.loc[:,"always_keep"] = False
    out = shrink_dataset(tmp_df,seqs_per_column=1/160)
    assert np.sum(out.keep) == 2

    # This will be only out.always_keep -- lowest we can go
    out = shrink_dataset(df,max_seq_number=1)
    assert np.sum(out.keep) == 4

    bad_max_seq_number = [-1,0,"test",int,14.2]
    for b in bad_max_seq_number:
        print(b)
        with pytest.raises(ValueError):
            shrink_dataset(df,max_seq_number=b)

    ## XX THIS IS A PARTIAL TEST AT THE MOMENT. FOCUSING ON TESTING LOWER LEVEL
    # FUNCTIONS
