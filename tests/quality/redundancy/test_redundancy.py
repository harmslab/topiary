
import pytest

import topiary
from topiary.quality import remove_redundancy, find_cutoff
from topiary.quality.redundancy.redundancy import _get_quality_scores, _reduce_redundancy_thread_manager
from topiary.quality.redundancy.redundancy import _EXPECTED_COLUMNS
from topiary.quality.redundancy.redundancy import _DummyTqdm

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
    assert _get_quality_scores(df.loc[0,:])[1] == 1
    assert _get_quality_scores(df.loc[0,:],{"Not a species":None})[1] == 1

    # Key species
    assert _get_quality_scores(df.loc[0,:],{species_in_df[0]:None})[1] == 0

    # Make sure quality assignment doing what we think
    for i in df.index:
        row = df.loc[i,:]
        scores = _get_quality_scores(row)

        values_from_df = list(row[_EXPECTED_COLUMNS])
        values_from_df.append(1/len(row.sequence))
        values_from_df = np.array(values_from_df,dtype=float)

        assert np.array_equal(scores[2:],values_from_df)

    # Test always_keep
    # No always_keep in datafframe
    assert _get_quality_scores(df.loc[0,:])[0] == 1

    df["always_keep"] = False # always keep False
    assert _get_quality_scores(df.loc[0,:])[0] == 1

    df["always_keep"] = True # always keep True
    assert _get_quality_scores(df.loc[0,:])[0] == 0

def test__reduce_redundancy_thread_manager():

    sequence_array = np.array(["STARE" for _ in range(4)])
    quality_array = np.array([np.zeros(3,dtype=int) for _ in range(4)])
    keep_array = np.ones(4,dtype=bool)
    cutoff = 0.9
    discard_key = True

    num_threads, all_args = _reduce_redundancy_thread_manager(sequence_array=sequence_array,
                                                              quality_array=quality_array,
                                                              keep_array=keep_array,
                                                              cutoff=cutoff,
                                                              discard_key=discard_key,
                                                              num_threads=1)

    assert num_threads == 1
    assert np.array_equal(all_args[0][0],(0,4))
    assert np.array_equal(all_args[0][1],(0,4))

    out = []
    for a in all_args:
        i_block = a[0]
        j_block = a[1]
        for i in range(i_block[0],i_block[1]):
            for j in range(j_block[0],j_block[1]):
                if i >= j:
                    continue
                out.append((i,j))

    out.sort()
    assert np.array_equal(out,((0,1),(0,2),(0,3),(1,2),(1,3),(2,3)))


    # Not worth chopping up problem for this small of an array -- should set
    # number of threads to 1.
    num_threads, all_args = _reduce_redundancy_thread_manager(sequence_array=sequence_array,
                                                              quality_array=quality_array,
                                                              keep_array=keep_array,
                                                              cutoff=cutoff,
                                                              discard_key=discard_key,
                                                              num_threads=2)

    assert num_threads == 1
    assert np.array_equal(all_args[0][0],(0,4))
    assert np.array_equal(all_args[0][1],(0,4))


    # assert np.array_equal(all_args[0][0],(0,2))
    # assert np.array_equal(all_args[0][1],(0,2))
    # assert np.array_equal(all_args[1][0],(0,2))
    # assert np.array_equal(all_args[1][1],(2,4))
    # assert np.array_equal(all_args[2][0],(2,4))
    # assert np.array_equal(all_args[2][1],(2,4))

    # sequence_array = np.array(["STARE" for _ in range(5)])
    # num_threads, all_args = _reduce_redundancy_thread_manager(sequence_array=sequence_array,
    #                                                           quality_array=quality_array,
    #                                                           keep_array=keep_array,
    #                                                           cutoff=cutoff,
    #                                                           discard_key=discard_key,
    #                                                           num_threads=20)
    #
    # assert num_threads == 1
    # print(all_args)
    # assert False
    # assert np.array_equal(all_args[0][0],(0,2))
    # assert np.array_equal(all_args[0][1],(0,2))
    # assert np.array_equal(all_args[1][0],(0,2))
    # assert np.array_equal(all_args[1][1],(2,4))
    # assert np.array_equal(all_args[2][0],(2,4))
    # assert np.array_equal(all_args[2][1],(2,4))


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

    bad_silent = [None,"test",int,float,{"test":1}]
    for b in bad_silent:
        print(f"trying bad silent {b}")
        with pytest.raises(ValueError):
            remove_redundancy(df=df,silent=b)

    good_silent = [True,False,0,1]
    for g in good_silent:
        print(f"trying good silent {g}")
        remove_redundancy(df=df,silent=g)


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
    assert np.sum(out_df.keep) < np.sum(df.keep)

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
    out_df = remove_redundancy(df=df,cutoff=0.2)
    assert np.sum(out_df.keep) == 1

def test_find_cutoff(test_dataframes):

    df = test_dataframes["good-df"].copy()

    # Should work
    find_cutoff(df=df)

    # -------------------------------------------------------------------------
    # Test argument parsing

    bad_df = [None,-1,1.1,"test",int,float,{"test":1},pd.DataFrame({"test":[1,2,3]})]
    for b in bad_df:
        with pytest.raises(ValueError):
            find_cutoff(df=b)

    # Bad min_cutoff
    bad_cutoff = [None,-1,1.1,"test",int,float,{"test":1},pd.DataFrame({"test":[1,2,3]})]
    for b in bad_cutoff:
        with pytest.raises(ValueError):
            find_cutoff(df=df,min_cutoff=b)

    # Bad max_cutoffs
    bad_cutoff = [None,-1,1.1,"test",int,float,{"test":1},pd.DataFrame({"test":[1,2,3]})]
    for b in bad_cutoff:
        with pytest.raises(ValueError):
            find_cutoff(df=df,max_cutoff=b)

    # Throw error because max smaller than min
    with pytest.raises(ValueError):
        find_cutoff(df=df,min_cutoff=0.9,max_cutoff=0.1)

    # Bad try_n_values
    bad_value = [None,-1,0,1.1,"test",int,float,{"test":1},pd.DataFrame({"test":[1,2,3]})]
    for b in bad_value:
        with pytest.raises(ValueError):
            find_cutoff(df=df,try_n_values=b)

    # Bad target_number
    bad_value = [None,-1,0,1.1,"test",int,float,{"test":1},pd.DataFrame({"test":[1,2,3]})]
    for b in bad_value:
        with pytest.raises(ValueError):
            find_cutoff(df=df,target_number=b)

    # Bad sample_size
    bad_value = [None,-1,0,1.1,"test",int,float,{"test":1},pd.DataFrame({"test":[1,2,3]})]
    for b in bad_value:
        with pytest.raises(ValueError):
            find_cutoff(df=df,sample_size=b)

    # Bad key species
    bad_key_species = [None,-1,1.1,"test",int,float,{"test":1}]
    for b in bad_key_species:
        print(f"trying bad key species {b}")
        with pytest.raises(ValueError):
            find_cutoff(df=df,key_species=b)

    # Make sure cutoffs work
    df = test_dataframes["good-df"].copy()
    for i in range(5):
        cutoff = find_cutoff(df,min_cutoff=0.5,max_cutoff=1.0,target_number=(i+1))
        new_df = remove_redundancy(df,cutoff=cutoff)

def test__DummyTqdm():

    pass

def test__DummyTqdm___enter__():

    pass

def test__FakeLock():

    pass
