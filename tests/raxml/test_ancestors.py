import pytest

import topiary
from topiary.raxml.ancestors import _make_ancestor_summary_trees
from topiary.raxml.ancestors import _parse_raxml_anc_output
from topiary.raxml.ancestors import generate_ancestors
from topiary.raxml.ancestors import _get_ancestral_gaps

import numpy as np

import os, ete3

def test__get_ancestral_gaps(tmpdir):

    # Check a variety of gapping
    alignment = \
    """6 9

    AAAAAAAAAA
    Q--QQ-Q-Q
    BBBBBBBBBB
    Q--QQ--Q-
    CCCCCCCCCC
    Q-Q-Q--QQ
    DDDDDDDDDD
    Q-Q-Q--Q-
    EEEEEEEEEE
    Q-QQ-Q-QQ
    FFFFFFFFFF
    Q-QQ-Q-Q-
    """

    os.mkdir(os.path.join(tmpdir,"toy-anc-data"))
    ali_file = os.path.join(tmpdir,"toy-anc-data","alignment.phy")
    tree_file = os.path.join(tmpdir,"toy-anc-data","tree.newick")

    f = open(ali_file,"w")
    f.write(alignment)
    f.close()

    tree = "(((AAAAAAAAAA,BBBBBBBBBB)AB,(CCCCCCCCCC,DDDDDDDDDD)CD)ABCD,(EEEEEEEEEE,FFFFFFFFFF)EF)ABCDEF;"
    T = ete3.Tree(tree,format=8)
    T.write(outfile=tree_file,format=8)

    gapping = _get_ancestral_gaps(ali_file,tree_file)

    expected = {'': [False, True, False, False, None, None, True, False, None],
                'ABCD': [False, True, False, False, False, True, True, False, None],
                'AB': [False, True, True, False, False, True, True, False, None],
                'CD': [False, True, False, True, False, True, True, False, None],
                'EF': [False, True, False, False, True, False, True, False, None]}

    for k in expected:
        assert np.array_equal(gapping[k],expected[k])

def test__make_ancestor_summary_trees():

    pass

def test__parse_raxml_anc_output():

    pass

def test_generate_ancestors():

    pass
