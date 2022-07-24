

import pytest

import topiary
from topiary.external.ncbi.blast.merge import _check_merge, merge_blast_df

import numpy as np
import pandas as pd

def test__check_merge():

    to_merge = [None,None,None]
    index = 1
    merge_index = 0
    expected = [0,0,None]
    to_merge = _check_merge(index,merge_index,to_merge)
    for i in range(len(to_merge)):
        assert to_merge[i] == expected[i]

    to_merge = [3,None,3]
    index = 1
    merge_index = 0

    expected = [3,3,3]
    to_merge = _check_merge(index,merge_index,to_merge)
    for i in range(len(to_merge)):
        assert to_merge[i] == expected[i]

    to_merge = [3,4,3,4]
    index = 1
    merge_index = 0

    expected = [3,3,3,3]
    to_merge = _check_merge(index,merge_index,to_merge)
    for i in range(len(to_merge)):
        assert to_merge[i] == expected[i]

    to_merge = [3,4,3,5,5]
    index = 1
    merge_index = 0

    expected = [3,3,3,5,5]
    to_merge = _check_merge(index,merge_index,to_merge)
    for i in range(len(to_merge)):
        assert to_merge[i] == expected[i]


def test_merge_blast_df(recip_blast_hit_dfs):

    blast_dfs = recip_blast_hit_dfs["ncbi"]

    concat = pd.concat(blast_dfs,ignore_index=True)
    merged = merge_blast_df(blast_dfs)

    # For these dataframes, all overlapping -- should merge perfectly
    assert len(np.unique(concat.accession)) == len(merged)

    df_0 = blast_dfs[0].copy()
    df_1 = blast_dfs[1].copy()

    merged = merge_blast_df([df_0,df_1])

    # Offset df_0 so the sequences no longer overlap. These should not merge.
    df_0.loc[:,"subject_start"] = df_0.subject_start + 1000
    df_0.loc[:,"subject_end"] = df_0.subject_end + 1000

    merged_offset = merge_blast_df([df_0,df_1])

    assert len(merged) == 15
    assert len(merged_offset) == len(df_0) + len(df_1)

    # Send in garbage
    bad_inputs = [["not","df"],[pd.DataFrame(),pd.DataFrame()],1,None,list,1.5]
    for b in bad_inputs:
        print(f"sending in bad input {b}")
        with pytest.raises(ValueError):
            merge_blast_df(b)

    # send in subtle garbage; missing requied column
    df_0 = blast_dfs[0].copy()
    df_1 = blast_dfs[1].copy()
    df_0 = df_0.drop(columns="accession")
    with pytest.raises(ValueError):
        merge_blast_df([df_0,df_1])
