
import pytest

import numpy as np

from topiary.quality.redundancy._block import _compare_seqs
from topiary.quality.redundancy.redundancy import _EXPECTED_COLUMNS

def test__compare_seqs(test_dataframes):

    A_seq = "TEST"
    B_seq = "TAST"

    # Identical quals
    A_qual = np.zeros(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    B_qual = np.zeros(len(_EXPECTED_COLUMNS) + 2,dtype=float)

    # Neither are always keep sequences
    A_qual[0] = 1
    B_qual[0] = 1

    # Neither are key sequences
    A_qual[1] = 1
    B_qual[1] = 1

    # Keep both; below cutoff
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.9)
    assert a1 is True
    assert a2 is True

    # Keep A arbitrarily
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is True
    assert a2 is False

    # Now make A_qual score worse than B, so keep B
    A_qual[2] = 1
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is False
    assert a2 is True

    # Not set up qual scores so neither are key_species, B has earlier better
    # score than A
    A_qual = np.ones(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    B_qual = np.ones(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    A_qual[-1] = 0
    B_qual[-2] = 0

    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is False
    assert a2 is True

    # both key species, A worse than B. No always_keep
    A_qual = np.zeros(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    B_qual = np.zeros(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    A_qual[0] = 1
    B_qual[0] = 1
    A_qual[2] = 1

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

    # both always keep, but  beter. Should keep both
    A_qual = np.zeros(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    B_qual = np.zeros(len(_EXPECTED_COLUMNS) + 2,dtype=float)
    A_qual[2] = 1
    a1, a2 = _compare_seqs(A_seq,B_seq,A_qual,B_qual,0.5)
    assert a1 is True
    assert a2 is True
