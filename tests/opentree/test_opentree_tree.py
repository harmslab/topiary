
import pytest

import topiary
from topiary.opentree.tree import get_species_tree
import numpy as np
import pandas as pd

import ete3

def test_get_species_tree(test_dataframes,df_with_species_not_resolvable):

    df = test_dataframes["good-df"]

    expected_num_leaves = len(np.unique(df.species))

    # make sure the uid, ott, and species loaded correctly
    T, dropped = get_species_tree(df)
    assert len(T.get_leaves()) == expected_num_leaves
    tips = [n.ott for n in T.get_leaves()]
    tips.sort()
    ott_from_df = list(df.loc[:,"ott"])
    ott_from_df.sort()
    assert np.array_equal(tips,ott_from_df)

    tips = [n.species for n in T.get_leaves()]
    tips.sort()
    from_df = list(df.loc[:,"species"])
    from_df.sort()
    assert np.array_equal(tips,from_df)

    # uid will be lists of length one since all species are unique in input
    # dataframe
    tips = [n.uid for n in T.get_leaves()]
    for t in tips:
        assert(len(t)) == 1
    tips = [t[0] for t in tips]
    tips.sort()
    from_df = list(df.loc[:,"uid"])
    from_df.sort()
    assert np.array_equal(tips,from_df)

    # nothing should have been dropped.
    assert len(dropped) == 0

    # Make sure the check for ott is working
    bad_df = df.drop(columns=["ott"])
    with pytest.raises(ValueError):
        T, dropped = get_species_tree(bad_df)

    # Test for resolvability check and for dropping keep = False
    T, dropped = get_species_tree(df_with_species_not_resolvable)

    assert len(dropped) == 1
    assert len(list(T.get_leaves())) == 357

    with pytest.raises(ValueError):
        T, dropped = get_species_tree(df_with_species_not_resolvable,
                                                            strict=True)
