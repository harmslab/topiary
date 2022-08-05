
import pytest

import topiary
from topiary.quality.alignment import _get_sparse_columns, _rle, _drop_gaps_only
#from topiary.quality.alignment import _find_too_many_sparse, _find_too_few_dense
#from topiary.quality.alignment import _find_long_insertions
from topiary.quality.alignment import AA_TO_INT, INT_TO_AA
from topiary.quality.alignment import score_alignment

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


def test_score_alignment():

    pass
