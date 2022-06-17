
import pytest

import topiary
from topiary.quality._taxonomic_sample import _prep_species_tree, _get_sequence_budgets

import ete3
import numpy as np
import pandas as pd

import random

def test__prep_species_tree(test_dataframes):

    df = test_dataframes["good-df"].copy()
    df["recip_paralog"] = "LY96"

    T = _prep_species_tree(df,paralog_column="recip_paralog")
    assert type(T) is ete3.Tree
    for leaf in T.get_leaves():
        assert np.array_equal(list(leaf.paralogs.keys()),["LY96"])
        assert len(leaf.paralogs["LY96"]) == 1
        assert leaf.uid[0] == leaf.paralogs["LY96"][0] # Make sure we're not mixing up uid

    with pytest.raises(ValueError):
        T = _prep_species_tree(df,paralog_column="not_a_column")

    df = test_dataframes["good-df"].copy()
    second_df = test_dataframes["good-df"].copy()

    df = pd.concat((df,second_df),ignore_index=True)
    df.uid = topiary._private.generate_uid(len(df))
    df["recip_paralog"] = "LY96"
    df.loc[0:4,"recip_paralog"] = "LY86"

    T = _prep_species_tree(df,paralog_column="recip_paralog")
    for leaf in T.get_leaves():
        assert np.array_equal(list(leaf.paralogs.keys()),["LY86","LY96"])
        assert len(leaf.paralogs["LY96"]) == 1
        assert len(leaf.paralogs["LY86"]) == 1


def test__get_sequence_budgets():

    some_tree = "(((A,B),(C,(D,H))),(E,F));"
    total_budget = 5
    CUTOFF = 1

    T = ete3.Tree(some_tree)

    for leaf in T.get_leaves():
        leaf.sequences = []
        if random.random() < CUTOFF:
            leaf.sequences.append("s")
            leaf.name = f"{leaf.name}*"

    for n in T.traverse():
        if hasattr(n,"budget"):
            n.__delattr__("budget")
        # if hasattr(n,"sequences"):
        #     for i in range(len(n.sequences)):
        #         n.sequences[i] = "s"
        #     if len(n.sequences) > 0:
        #         n.name = f"{n.name.split('_')[0]}_s"
        #     else:
        #         n.name = n.name.split('_')[0]


    pass
