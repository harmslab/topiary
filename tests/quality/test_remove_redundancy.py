
import pytest

import topiary
from topiary.quality import remove_redundancy
from topiary.quality.remove_redundancy import _get_quality_scores, _compare_seqs
from topiary.quality.remove_redundancy import _EXPECTED_COLUMNS, _LENGTH_COLUMN

import numpy as np
import pandas as pd

def test__get_quality_scores(test_dataframes):

    # Get copy of the dataframe -- we're going to hack it
    df = test_dataframes["good-df"].copy()
    species_in_df = list(df.species)

    # Needs to be defined.
    df["diff_from_median"] = 0

    # Make sure that the key species encoding works as expected
    # No key_species passed in or key_species not in dataframe
    assert _get_quality_scores(df.loc[0,:])[0] == 1
    assert _get_quality_scores(df.loc[0,:],{"Not a species":None})[0] == 1

    # Key species
    assert _get_quality_scores(df.loc[0,:],{species_in_df[0]:None})[0] == 0

    # Make sure quality assignment doing what we think
    for i in df.index:
        row = df.loc[i,:]
        scores = _get_quality_scores(row)

        values_from_df = np.array(row[_EXPECTED_COLUMNS],dtype=float)
        values_from_df[_LENGTH_COLUMN] = 1/values_from_df[_LENGTH_COLUMN]

        assert np.array_equal(scores[1:],values_from_df)

def test__compare_seqs(test_dataframes):

    A_seq = "TEST"
    B_seq = "TAST"

    # Identical quals
    A_qual = np.zeros(len(_EXPECTED_COLUMNS) + 1,dtype=float)
    B_qual = np.zeros(len(_EXPECTED_COLUMNS) + 1,dtype=float)

    # Neither are key sequences
    A_qual[0] = 1
    B_qual[0] = 1

    # Keep both; below cutoff
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.9)
    assert a1 is True
    assert a2 is True

    # Keep A arbitrarily
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is True
    assert a2 is False

    # Now make A_qual score worse than B, so keep B
    A_qual[1] = 1
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is False
    assert a2 is True

    # Not set up qual scores so neither are key_species, B has earlier better
    # score than A
    A_qual = np.ones(len(_EXPECTED_COLUMNS) + 1,dtype=float)
    B_qual = np.ones(len(_EXPECTED_COLUMNS) + 1,dtype=float)
    A_qual[-1] = 0
    B_qual[-2] = 0

    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is False
    assert a2 is True

    # both key species, A worse than B
    A_qual = np.zeros(len(_EXPECTED_COLUMNS) + 1,dtype=float)
    B_qual = np.zeros(len(_EXPECTED_COLUMNS) + 1,dtype=float)
    A_qual[1] = 1

    # implicit discard_key flag
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is True
    assert a2 is True

    # Explicit discard_key flag
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5,discard_key=False)
    assert a1 is True
    assert a2 is True

    # Check discard_key flag
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5,discard_key=True)
    assert a1 is False
    assert a2 is True

def test_remove_redundancy(test_dataframes):

    df = test_dataframes["good-df"].copy()

    # -------------------------------------------------------------------------
    # Test argument parsing

    bad_df = [None,-1,1.1,"test",int,float,{"test":1},pd.DataFrame({"test":[1,2,3]})]
    for b in bad_df:
        with pytest.raises(ValueError):
            remove_redundancy(df=b)

    remove_redundancy(df=df)

    bad_cutoff = [None,-1,1.1,"test",int,float,{"test":1},pd.DataFrame({"test":[1,2,3]})]
    for b in bad_cutoff:
        with pytest.raises(ValueError):
            remove_redundancy(df=df,cutoff=b)

    good_cutoff = [0,0.5,1]
    for g in good_cutoff:
        remove_redundancy(df=df,cutoff=g)

    bad_key_species = [None,-1,1.1,"test",int,float,{"test":1}]
    for b in bad_key_species:
        print(f"trying bad key species {b}")
        with pytest.raises(ValueError):
            remove_redundancy(df=df,key_species=b)

    good_key_species = [[],["test"],("test","this"),np.array(["test","this"])]
    for g in good_key_species:
        print(f"trying good key species {g}")
        remove_redundancy(df=df,key_species=g)


    # -------------------------------------------------------------------------
    # Make sure dropping is happening a sane way that depends on cutoff and
    # key_species.

    df = test_dataframes["good-df"].copy()
    species_in_df = list(df.species)

    # sequences in this dataframe are between 0.9125 and 0.98125 identical.
    out_df = remove_redundancy(df=df,cutoff=0.99)
    assert np.sum(out_df.keep) == np.sum(df.keep)

    # Cut some
    out_df = remove_redundancy(df=df,cutoff=0.96)
    assert np.sum(out_df.keep) <= np.sum(df.keep)

    # Cut basically all -- only one shoudl survive
    out_df = remove_redundancy(df=df,cutoff=0.50)
    assert np.sum(out_df.keep) == 1

    # All key species -- all should survive
    out_df = remove_redundancy(df=df,cutoff=0.50,key_species=species_in_df)
    assert np.sum(out_df.keep) == np.sum(df.keep)

    # One isn't keep -- make sure it's dropped
    out_df = remove_redundancy(df=df,cutoff=0.50,key_species=species_in_df[1:])
    assert np.sum(out_df.keep) == 4
    assert out_df.loc[out_df["species"] == species_in_df[0],:].iloc[0].keep == False

    # Now make all the same species -- should get all dropped even if key_species
    df.species = species_in_df[0]
    out_df = remove_redundancy(df=df,cutoff=0.50,key_species=[species_in_df[0]])
    assert np.sum(out_df.keep) == 1

    # shouldn't matter what is in key_species here
    df.species = species_in_df[0]
    out_df = remove_redundancy(df=df,cutoff=0.50,key_species=[])
    assert np.sum(out_df.keep) == 1

    # -------------------------------------------------------------------------
    # Make sure it takes row with higher quality

    df = test_dataframes["good-df"].copy()
    df = df.iloc[2:4]
    out_df = remove_redundancy(df=df,cutoff=0.50)
    assert out_df.keep.iloc[0] == False
    assert out_df.keep.iloc[1] == True

    df = test_dataframes["good-df"].copy()
    df = df.iloc[1:3]
    out_df = remove_redundancy(df=df,cutoff=0.50)
    assert out_df.keep.iloc[0] == True
    assert out_df.keep.iloc[1] == False

    # Make sure length bit is being processed properly. Make first sequence short
    # so it gets dropped
    df = test_dataframes["good-df_only-required-columns"].copy()
    df.loc[0,"sequence"] = "MLPFLFFS"
    df.loc[:,"length"] = [len(s) for s in df.loc[:,"sequence"]]
    out_df = remove_redundancy(df=df,cutoff=0.50)
    assert out_df.keep.iloc[0] == False

    # -------------------------------------------------------------------------
    # Check input dataframe without quality information

    # Should work fine but cut nothing
    df = test_dataframes["good-df_only-required-columns"].copy()
    out_df = remove_redundancy(df=df,cutoff=0.99)
    assert np.sum(out_df.keep) == len(df.sequence)

    # Cut some
    out_df = remove_redundancy(df=df,cutoff=0.96)
    assert np.sum(out_df.keep) <= len(df.sequence)

    # Cut basically all -- only one should survive
    out_df = remove_redundancy(df=df,cutoff=0.50)
    assert np.sum(out_df.keep) == 1
