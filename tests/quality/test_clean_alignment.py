
import pytest

import topiary
from topiary.quality._clean_alignment import _get_sparse_columns, _rle, _drop_gaps_only
from topiary.quality._clean_alignment import _find_too_many_sparse, _find_too_few_dense
from topiary.quality._clean_alignment import _find_long_insertions, clean_alignment
from topiary.quality._clean_alignment import AA_TO_INT, INT_TO_AA

import numpy as np
import pandas as pd

def data_for_test():

    good =  "TEST---------------S-"
    long =  "TESTTHISTHISTHISTHIS-"
    short = "------------------IS-"

    seqs = []
    for i in range(20):
        seqs.append([AA_TO_INT[g] for g in list(good)])
    seqs.append([AA_TO_INT[g] for g in list(long)])
    seqs.append([AA_TO_INT[g] for g in list(short)])

    return np.array(seqs)

def test__get_sparse_columns():

    seqs = data_for_test()

    sparse_mask = np.ones(seqs.shape[1],dtype=bool)
    sparse_mask[:4] = 0
    sparse_mask[-2] = 0

    sparse = _get_sparse_columns(seqs)
    assert np.array_equal(sparse,sparse_mask)

    sparse = _get_sparse_columns(seqs,sparse_column_cutoff=0.045)
    sparse_mask = np.ones(seqs.shape[1],dtype=bool)
    sparse_mask[-2] = 0
    assert np.array_equal(sparse,sparse_mask)

def test__rle():

    input = np.array([1,1,1,0,0,0,1,1,1])
    run_lengths, start_positions, values = _rle(input)
    assert np.array_equal(run_lengths,np.array([3,3,3]))
    assert np.array_equal(start_positions,np.array([0,3,6]))
    assert np.array_equal(values,np.array([1,0,1]))

    input = np.array([1,2,1,0])
    run_lengths, start_positions, values = _rle(input)
    assert np.array_equal(run_lengths,np.array([1,1,1,1]))
    assert np.array_equal(start_positions,np.array([0,1,2,3]))
    assert np.array_equal(values,np.array([1,2,1,0]))

    input = np.array([0,0,0,0])
    run_lengths, start_positions, values = _rle(input)
    assert np.array_equal(run_lengths,np.array([4]))
    assert np.array_equal(start_positions,np.array([0]))
    assert np.array_equal(values,np.array([0]))

def test__drop_gaps_only():

    seqs = data_for_test()

    new_seqs = _drop_gaps_only(seqs)
    assert new_seqs.shape[1] == seqs.shape[1] - 1

    hacked_seqs = seqs[:20]
    new_seqs = _drop_gaps_only(hacked_seqs)
    assert new_seqs.shape[1] == 5

def test__find_too_many_sparse():

    seqs = data_for_test()
    cumulative_keep = np.arange(seqs.shape[0],dtype=int)
    force_keep = np.zeros(seqs.shape[0],dtype=bool)

    # Should drop one long insertion
    new_seqs, new_cum_keep, new_force = _find_too_many_sparse(seqs,
                                                              cumulative_keep,
                                                              force_keep=force_keep,
                                                              sparse_column_cutoff=0.9,
                                                              maximum_sparse_allowed=0.5)
    assert new_seqs.shape[0] == seqs.shape[0] - 1
    assert new_seqs.shape[1] == 6
    new_keep_mask = np.ones(seqs.shape[0],dtype=bool)
    new_keep_mask[-2] = False
    assert np.array_equal(new_cum_keep,cumulative_keep[new_keep_mask])

    # Should not drop -- requires too much of sequence be in insertion
    new_seqs, new_cum_keep, new_force = _find_too_many_sparse(seqs,
                                                              cumulative_keep,
                                                              force_keep=force_keep,
                                                              sparse_column_cutoff=0.9,
                                                              maximum_sparse_allowed=0.9)
    assert new_seqs.shape[0] == seqs.shape[0]
    assert new_seqs.shape[1] == seqs.shape[1] - 1 # drop gap at end
    new_keep_mask = np.ones(seqs.shape[0],dtype=bool)
    assert np.array_equal(new_cum_keep,cumulative_keep[new_keep_mask])

    # Should not drop -- requires very strict sparse sequence call
    new_seqs, new_cum_keep, new_force = _find_too_many_sparse(seqs,
                                                              cumulative_keep,
                                                              force_keep=force_keep,
                                                              sparse_column_cutoff=0.96,
                                                              maximum_sparse_allowed=0.5)

    assert new_seqs.shape[0] == seqs.shape[0]
    assert new_seqs.shape[1] == seqs.shape[1] - 1 # drop gap at end
    new_keep_mask = np.ones(seqs.shape[0],dtype=bool)
    assert np.array_equal(new_cum_keep,cumulative_keep[new_keep_mask])

    # Should not drop -- force_keep set to True for the sequence we would drop
    force_keep[-2] = True
    new_seqs, new_cum_keep, new_force = _find_too_many_sparse(seqs,
                                                              cumulative_keep,
                                                              force_keep=force_keep,
                                                              sparse_column_cutoff=0.9,
                                                              maximum_sparse_allowed=0.5)

    assert new_seqs.shape[0] == seqs.shape[0]
    assert new_seqs.shape[1] == seqs.shape[1] - 1 # drop gap at end
    new_keep_mask = np.ones(seqs.shape[0],dtype=bool)
    assert np.array_equal(new_cum_keep,cumulative_keep[new_keep_mask])

def test__find_too_few_dense():

    seqs = data_for_test()
    cumulative_keep = np.arange(seqs.shape[0],dtype=int)
    force_keep = np.zeros(seqs.shape[0],dtype=bool)

    # Should drop sequence that is super short
    new_seqs, new_cum_keep, new_force = _find_too_few_dense(seqs,
                                                            cumulative_keep,
                                                            force_keep=force_keep,
                                                            sparse_column_cutoff=0.9,
                                                            minimum_dense_required=0.9)
    assert new_seqs.shape[0] == seqs.shape[0] - 1
    assert new_seqs.shape[1] == seqs.shape[1] - 1 # drop gap at end
    new_keep_mask = np.ones(seqs.shape[0],dtype=bool)
    new_keep_mask[-1] = False
    assert np.array_equal(new_cum_keep,cumulative_keep[new_keep_mask])

    # Should not drop short sequence -- has enough of dense
    new_seqs, new_cum_keep, new_force = _find_too_few_dense(seqs,
                                                            cumulative_keep,
                                                            force_keep=force_keep,
                                                            sparse_column_cutoff=0.9,
                                                            minimum_dense_required=0.1)
    assert new_seqs.shape[0] == seqs.shape[0]
    assert new_seqs.shape[1] == seqs.shape[1] - 1 # drop gap at end
    new_keep_mask = np.ones(seqs.shape[0],dtype=bool)
    assert np.array_equal(new_cum_keep,cumulative_keep[new_keep_mask])

    # Should drop every sequence but super long sequence. Every column called
    # as dense, so every other sequence does not have enough fx dense
    new_seqs, new_cum_keep, new_force = _find_too_few_dense(seqs,
                                                            cumulative_keep,
                                                            force_keep=force_keep,
                                                            sparse_column_cutoff=0.96,
                                                            minimum_dense_required=0.9)
    assert new_seqs.shape[0] == 1
    assert new_seqs.shape[1] == seqs.shape[1] - 1 # drop gap at end
    new_keep_mask = np.zeros(seqs.shape[0],dtype=bool)
    new_keep_mask[-2] = True
    assert np.array_equal(new_cum_keep,cumulative_keep[new_keep_mask])

    # Should not drop any sequence. Set force_keep to True for all sequences
    force_keep = np.ones(seqs.shape[0],dtype=bool)
    new_seqs, new_cum_keep, new_force = _find_too_few_dense(seqs,
                                                            cumulative_keep,
                                                            force_keep=force_keep,
                                                            sparse_column_cutoff=0.96,
                                                            minimum_dense_required=0.9)
    assert new_seqs.shape[0] == seqs.shape[0]
    assert new_seqs.shape[1] == seqs.shape[1] - 1 # drop gap at end
    new_keep_mask = np.ones(seqs.shape[0],dtype=bool)
    assert np.array_equal(new_cum_keep,cumulative_keep[new_keep_mask])

def test__find_long_insertions():

    seqs = data_for_test()
    cumulative_keep = np.arange(seqs.shape[0],dtype=int)
    force_keep = np.zeros(seqs.shape[0],dtype=bool)

    # Should drop long sequence
    new_seqs, new_cum_keep, new_force = _find_long_insertions(seqs,
                                                              cumulative_keep,
                                                              force_keep=force_keep,
                                                              sparse_column_cutoff=0.9,
                                                              long_insertion_length=5,
                                                              long_insertion_fx_cutoff=0.9)
    assert new_seqs.shape[0] == seqs.shape[0] - 1
    assert new_seqs.shape[1] == 6 # Drop massive gaps because we dropped long seq
    new_keep_mask = np.ones(seqs.shape[0],dtype=bool)
    new_keep_mask[-2] = False
    assert np.array_equal(new_cum_keep,cumulative_keep[new_keep_mask])

    # Should not drop long sequence -- looking for gap that is too long
    new_seqs, new_cum_keep, new_force = _find_long_insertions(seqs,
                                                              cumulative_keep,
                                                              force_keep=force_keep,
                                                              sparse_column_cutoff=0.9,
                                                              long_insertion_length=16,
                                                              long_insertion_fx_cutoff=0.9)
    assert new_seqs.shape[0] == seqs.shape[0]
    assert new_seqs.shape[1] == seqs.shape[1] - 1 # drop gap at end
    new_keep_mask = np.ones(seqs.shape[0],dtype=bool)
    assert np.array_equal(new_cum_keep,cumulative_keep[new_keep_mask])

    # Should not drop long sequence -- looking for too much gap fx. (This is
    # silly because fx_cutoff will always be 1.0 or less when really run;
    # for testing purposes, drop in something wacky)
    new_seqs, new_cum_keep, new_force = _find_long_insertions(seqs,
                                                              cumulative_keep,
                                                              force_keep=force_keep,
                                                              sparse_column_cutoff=0.9,
                                                              long_insertion_length=5,
                                                              long_insertion_fx_cutoff=1.1)
    assert new_seqs.shape[0] == seqs.shape[0]
    assert new_seqs.shape[1] == seqs.shape[1] - 1 # drop gap at end
    new_keep_mask = np.ones(seqs.shape[0],dtype=bool)
    assert np.array_equal(new_cum_keep,cumulative_keep[new_keep_mask])

    # Should not drop long sequence -- set force_keep to True for that sequence
    force_keep[-2] = True
    new_seqs, new_cum_keep, new_force = _find_long_insertions(seqs,
                                                              cumulative_keep,
                                                              force_keep=force_keep,
                                                              sparse_column_cutoff=0.9,
                                                              long_insertion_length=5,
                                                              long_insertion_fx_cutoff=1.1)
    assert new_seqs.shape[0] == seqs.shape[0]
    assert new_seqs.shape[1] == seqs.shape[1] - 1 # drop gap at end
    new_keep_mask = np.ones(seqs.shape[0],dtype=bool)
    assert np.array_equal(new_cum_keep,cumulative_keep[new_keep_mask])

def test_clean_alignment(test_dataframes):

    pass
