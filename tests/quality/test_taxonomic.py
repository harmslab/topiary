
import pytest

import topiary
from topiary.quality import taxonomic as tx

import ete3
import numpy as np
import pandas as pd

import random

def test__prep_species_tree(test_dataframes):

    df = test_dataframes["good-df"].copy()
    df["recip_paralog"] = "LY96"

    T = tx._prep_species_tree(df,paralog_column="recip_paralog")
    assert type(T) is ete3.Tree
    for leaf in T.get_leaves():
        assert np.array_equal(list(leaf.paralogs.keys()),["LY96"])
        assert len(leaf.paralogs["LY96"]) == 1
        assert leaf.uid[0] == leaf.paralogs["LY96"][0] # Make sure we're not mixing up uid

    with pytest.raises(ValueError):
        T = tx._prep_species_tree(df,paralog_column="not_a_column")

    df = test_dataframes["good-df"].copy()
    second_df = test_dataframes["good-df"].copy()

    df = pd.concat((df,second_df),ignore_index=True)
    df.uid = topiary._private.generate_uid(len(df))
    df["recip_paralog"] = "LY96"
    df.loc[0:4,"recip_paralog"] = "LY86"

    T = tx._prep_species_tree(df,paralog_column="recip_paralog")
    for leaf in T.get_leaves():
        assert np.array_equal(list(leaf.paralogs.keys()),["LY86","LY96"])
        assert len(leaf.paralogs["LY96"]) == 1
        assert len(leaf.paralogs["LY86"]) == 1

def test__even_paralog_budgeting():

    some_tree = "(((A,B),(C,(D,H))),(E,F));"

    # Make tree with X and Y paralogs at the tip, each of which are seen once
    T = ete3.Tree(some_tree)
    for leaf in T.get_leaves():
        leaf.paralogs = {}
        leaf.paralogs["X"] = [topiary._private.generate_uid(1)]
        leaf.paralogs["Y"] = [topiary._private.generate_uid(1)]
    T_bak = T.copy()

    T = T_bak.copy()
    budget = tx._even_paralog_budgeting(T,overall_budget=18)
    assert budget["X"] == 7
    assert budget["Y"] == 7

    T = T_bak.copy()
    budget = tx._even_paralog_budgeting(T,overall_budget=10)
    assert budget["X"] == 5
    assert budget["Y"] == 5

    T = T_bak.copy()
    budget = tx._even_paralog_budgeting(T,overall_budget=0)
    assert budget["X"] == 0
    assert budget["Y"] == 0

    T = T_bak.copy()
    budget = tx._even_paralog_budgeting(T,overall_budget=1000)
    assert budget["X"] == 7
    assert budget["Y"] == 7

    T = ete3.Tree(some_tree)
    for leaf in T.get_leaves():
        leaf.paralogs = {}
        leaf.paralogs["X"] = [topiary._private.generate_uid(1)]
        leaf.paralogs["Y"] = topiary._private.generate_uid(2)
    T_bak = T.copy()

    T = T_bak.copy()
    budget = tx._even_paralog_budgeting(T,overall_budget=21)
    assert budget["X"] == 7
    assert budget["Y"] == 14

    T = T_bak.copy()
    budget = tx._even_paralog_budgeting(T,overall_budget=10)
    assert budget["X"] == 5
    assert budget["Y"] == 5

    T = T_bak.copy()
    budget = tx._even_paralog_budgeting(T,overall_budget=4)
    assert budget["X"] == 2
    assert budget["Y"] == 2

    T = T_bak.copy()
    budget = tx._even_paralog_budgeting(T,overall_budget=500)
    assert budget["X"] == 7
    assert budget["Y"] == 14

    T = T_bak.copy()
    budget = tx._even_paralog_budgeting(T,overall_budget=0)
    assert budget["X"] == 0
    assert budget["Y"] == 0


def test__weighted_paralog_budgeting():

    some_tree = "(((A,B),(C,(D,H))),(E,F));"

    # Make tree with X and Y paralogs at the tip, each of which are seen once
    T = ete3.Tree(some_tree)
    for leaf in T.get_leaves():
        leaf.paralogs = {}
        leaf.paralogs["X"] = [topiary._private.generate_uid(1)]
        leaf.paralogs["Y"] = [topiary._private.generate_uid(1)]
    T_bak = T.copy()

    T = T_bak.copy()
    budget = tx._weighted_paralog_budgeting(T,overall_budget=18)
    assert budget["X"] == 7
    assert budget["Y"] == 7

    T = T_bak.copy()
    budget = tx._weighted_paralog_budgeting(T,overall_budget=10)
    assert budget["X"] == 5
    assert budget["Y"] == 5

    T = T_bak.copy()
    budget = tx._weighted_paralog_budgeting(T,overall_budget=0)
    assert budget["X"] == 0
    assert budget["Y"] == 0

    T = T_bak.copy()
    budget = tx._weighted_paralog_budgeting(T,overall_budget=1000)
    assert budget["X"] == 7
    assert budget["Y"] == 7

    T = ete3.Tree(some_tree)
    for leaf in T.get_leaves():
        leaf.paralogs = {}
        leaf.paralogs["X"] = [topiary._private.generate_uid(1)]
        leaf.paralogs["Y"] = topiary._private.generate_uid(2)
    T_bak = T.copy()

    T = T_bak.copy()
    budget = tx._weighted_paralog_budgeting(T,overall_budget=21)
    assert budget["X"] == 7
    assert budget["Y"] == 14

    T = T_bak.copy()
    budget = tx._weighted_paralog_budgeting(T,overall_budget=9)
    assert budget["X"] == 3
    assert budget["Y"] == 6

    T = T_bak.copy()
    budget = tx._weighted_paralog_budgeting(T,overall_budget=4)
    assert budget["X"] == 1
    assert budget["Y"] == 3

    T = T_bak.copy()
    budget = tx._weighted_paralog_budgeting(T,overall_budget=500)
    assert budget["X"] == 7
    assert budget["Y"] == 14

    T = T_bak.copy()
    budget = tx._weighted_paralog_budgeting(T,overall_budget=0)
    assert budget["X"] == 0
    assert budget["Y"] == 0

def test__finalize_paralog_budget():

    paralog_budget = {'Y': 10, 'X': 11}
    paralog_counts = {'Y': 14, 'X': 7}
    final_budget = tx._finalize_paralog_budget(paralog_budget,paralog_counts)
    assert final_budget["X"] == 7
    assert final_budget["Y"] == 14

    paralog_budget = {'Y': 10, 'X': 11}
    paralog_counts = {'Y': 7, 'X': 7}
    final_budget = tx._finalize_paralog_budget(paralog_budget,paralog_counts)
    assert final_budget["X"] == 7
    assert final_budget["Y"] == 7

    paralog_budget = {'Y': 10, 'X': 0}
    paralog_counts = {'Y': 0, 'X': 10}
    final_budget = tx._finalize_paralog_budget(paralog_budget,paralog_counts)
    assert final_budget["X"] == 10
    assert final_budget["Y"] == 0

    paralog_budget = {'Y': 5, 'X': 5}
    paralog_counts = {'Y': 6, 'X': 5}
    final_budget = tx._finalize_paralog_budget(paralog_budget,paralog_counts)
    assert final_budget["X"] == 5
    assert final_budget["Y"] == 5

    paralog_budget = {'Y': 5, 'X': 1}
    paralog_counts = {'Y': 6, 'X': 10000}
    final_budget = tx._finalize_paralog_budget(paralog_budget,paralog_counts)
    assert final_budget["X"] == 1
    assert final_budget["Y"] == 5

def test__get_sequence_budgets():

    some_tree = "(((A,B),(C,(D,H))),(E,F));"

    # Make tree with single sequence at each tip
    T = ete3.Tree(some_tree)
    for leaf in T.get_leaves():
        leaf.sequences = ("X",)
    T_bak = T.copy()

    # Keep every tip
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,7)
    for leaf in T.get_leaves():
        assert leaf.budget == 1
    assert T.get_common_ancestor(("A","B")).budget == 2
    assert T.get_common_ancestor(("D","H")).budget == 2
    assert T.get_common_ancestor(("D","H","C")).budget == 3
    assert T.get_common_ancestor(("E","F")).budget == 2
    assert T.get_common_ancestor(("D","H","C","B","A")).budget == 5
    assert T.get_common_ancestor(("D","H","C","B","A","E","F")).budget == 7

    # Keep subset of tips
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,5)
    expected = {"A":0,"B":0,"C":1,"D":0,"H":0,"E":1,"F":1}
    for leaf in T.get_leaves():
        assert leaf.budget == expected[leaf.name]
    assert T.get_common_ancestor(("A","B")).budget == 1
    assert T.get_common_ancestor(("D","H")).budget == 1
    assert T.get_common_ancestor(("D","H","C")).budget == 2
    assert T.get_common_ancestor(("E","F")).budget == 2
    assert T.get_common_ancestor(("D","H","C","B","A")).budget == 3
    assert T.get_common_ancestor(("D","H","C","B","A","E","F")).budget == 5

    # Keep smaller subset of tips
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,4)
    expected = {"A":0,"B":0,"C":0,"D":0,"H":0,"E":1,"F":1}
    for leaf in T.get_leaves():
        assert leaf.budget == expected[leaf.name]
    assert T.get_common_ancestor(("A","B")).budget == 1
    assert T.get_common_ancestor(("D","H")).budget == 0
    assert T.get_common_ancestor(("D","H","C")).budget == 1
    assert T.get_common_ancestor(("E","F")).budget == 2
    assert T.get_common_ancestor(("D","H","C","B","A")).budget == 2
    assert T.get_common_ancestor(("D","H","C","B","A","E","F")).budget == 4

    # Keep smaller subset of tips
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,2)
    expected = {"A":0,"B":0,"C":0,"D":0,"H":0,"E":0,"F":0}
    for leaf in T.get_leaves():
        assert leaf.budget == expected[leaf.name]
    assert T.get_common_ancestor(("A","B")).budget == 0
    assert T.get_common_ancestor(("D","H")).budget == 0
    assert T.get_common_ancestor(("D","H","C")).budget == 0
    assert T.get_common_ancestor(("E","F")).budget == 1
    assert T.get_common_ancestor(("D","H","C","B","A")).budget == 1
    assert T.get_common_ancestor(("D","H","C","B","A","E","F")).budget == 2

    # Extra budget
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,500)
    for leaf in T.get_leaves():
        assert leaf.budget == 1
    assert T.get_common_ancestor(("A","B")).budget == 2
    assert T.get_common_ancestor(("D","H")).budget == 2
    assert T.get_common_ancestor(("D","H","C")).budget == 3
    assert T.get_common_ancestor(("E","F")).budget == 2
    assert T.get_common_ancestor(("D","H","C","B","A")).budget == 5
    assert T.get_common_ancestor(("D","H","C","B","A","E","F")).budget == 7

    # No budget
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,0)
    for leaf in T.get_leaves():
        assert leaf.budget == 0
    assert T.get_common_ancestor(("A","B")).budget == 0
    assert T.get_common_ancestor(("D","H")).budget == 0
    assert T.get_common_ancestor(("D","H","C")).budget == 0
    assert T.get_common_ancestor(("E","F")).budget == 0
    assert T.get_common_ancestor(("D","H","C","B","A")).budget == 0
    assert T.get_common_ancestor(("D","H","C","B","A","E","F")).budget == 0

    # Make tree with one tip that has multiple sequences
    some_tree = "(((A,B),(C,(D,H))),(E,F));"
    T = ete3.Tree(some_tree)
    for leaf in T.get_leaves():
        if leaf.name == "E":
            leaf.sequences = ("X","Y","Z")
        else:
            leaf.sequences = ("X",)
    T_bak = T.copy()

    # Keep every tip
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,9)
    expected = {"A":1,"B":1,"C":1,"D":1,"H":1,"E":3,"F":1}
    for leaf in T.get_leaves():
        assert leaf.budget == expected[leaf.name]
    assert T.get_common_ancestor(("A","B")).budget == 2
    assert T.get_common_ancestor(("D","H")).budget == 2
    assert T.get_common_ancestor(("D","H","C")).budget == 3
    assert T.get_common_ancestor(("E","F")).budget == 4
    assert T.get_common_ancestor(("D","H","C","B","A")).budget == 5
    assert T.get_common_ancestor(("D","H","C","B","A","E","F")).budget == 9

    # Keep subset of tips
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,5)
    expected = {"A":0,"B":0,"C":1,"D":0,"H":0,"E":1,"F":1}
    for leaf in T.get_leaves():
        assert leaf.budget == expected[leaf.name]
    assert T.get_common_ancestor(("A","B")).budget == 1
    assert T.get_common_ancestor(("D","H")).budget == 1
    assert T.get_common_ancestor(("D","H","C")).budget == 2
    assert T.get_common_ancestor(("E","F")).budget == 2
    assert T.get_common_ancestor(("D","H","C","B","A")).budget == 3
    assert T.get_common_ancestor(("D","H","C","B","A","E","F")).budget == 5

    # Keep subset of tips
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,3)
    expected = {"A":0,"B":0,"C":0,"D":0,"H":0,"E":0,"F":0}
    for leaf in T.get_leaves():
        assert leaf.budget == expected[leaf.name]
    assert T.get_common_ancestor(("A","B")).budget == 1
    assert T.get_common_ancestor(("D","H")).budget == 0
    assert T.get_common_ancestor(("D","H","C")).budget == 1
    assert T.get_common_ancestor(("E","F")).budget == 1
    assert T.get_common_ancestor(("D","H","C","B","A")).budget == 2
    assert T.get_common_ancestor(("D","H","C","B","A","E","F")).budget == 3

def test__identify_merge_blocks():

    some_tree = "(((A,B),(C,(D,H))),(E,F));"

    # Make a tree with sequences at the tips
    T = ete3.Tree(some_tree)
    for leaf in T.get_leaves():
        leaf.sequences = (topiary._private.uid.generate_uid(1),)
    T_bak = T.copy()


    # Keep everyone
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,7)
    merge_blocks = tx._identify_merge_blocks(T)
    expected_merges = {("A",):1,("B",):1,("C",):1,("D",):1,("H",):1,("E",):1,("F",):1}
    for m in merge_blocks:
        leaves = [leaf.name for leaf in m[2].get_leaves()]
        leaves.sort()
        leaves = tuple(leaves)
        assert m[0] == expected_merges[leaves]

    # Only keep 5
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,5)
    merge_blocks = tx._identify_merge_blocks(T)
    expected_merges = {("A","B"):1,("C",):1,("D","H",):1,("E",):1,("F",):1}
    for m in merge_blocks:
        leaves = [leaf.name for leaf in m[2].get_leaves()]
        leaves.sort()
        leaves = tuple(leaves)
        assert m[0] == expected_merges[leaves]

    # Only keep 3
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,3)
    merge_blocks = tx._identify_merge_blocks(T)
    expected_merges = {("A","B"):1,("C","D","H",):1,("E","F",):1}
    for m in merge_blocks:
        leaves = [leaf.name for leaf in m[2].get_leaves()]
        leaves.sort()
        leaves = tuple(leaves)
        assert m[0] == expected_merges[leaves]

    # Only keep 1
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,1)
    merge_blocks = tx._identify_merge_blocks(T)
    expected_merges = {("A","B","C","D","E","F","H",):1}
    for m in merge_blocks:
        leaves = [leaf.name for leaf in m[2].get_leaves()]
        leaves.sort()
        leaves = tuple(leaves)
        assert m[0] == expected_merges[leaves]


    some_tree = "(((A,B),(C,(D,H))),(E,F));"

    # Make a tree with sequences at the tips. Give E two sequences
    T = ete3.Tree(some_tree)
    for leaf in T.get_leaves():
        if leaf.name == "E":
            leaf.sequences = tuple(topiary._private.uid.generate_uid(2))
        else:
            leaf.sequences = (topiary._private.uid.generate_uid(1),)
    T_bak = T.copy()

    # Keep everyone
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,8)
    merge_blocks = tx._identify_merge_blocks(T)
    expected_merges = {("A",):1,("B",):1,("C",):1,("D",):1,("H",):1,("E",):2,("F",):1}
    for m in merge_blocks:
        leaves = [leaf.name for leaf in m[2].get_leaves()]
        leaves.sort()
        leaves = tuple(leaves)
        assert m[0] == expected_merges[leaves]

    # Keep everyone but one
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,7)
    merge_blocks = tx._identify_merge_blocks(T)
    expected_merges = {("A",):1,("B",):1,("C",):1,("D","H",):1,("E",):2,("F",):1}
    for m in merge_blocks:
        leaves = [leaf.name for leaf in m[2].get_leaves()]
        leaves.sort()
        leaves = tuple(leaves)
        assert m[0] == expected_merges[leaves]
        if leaves == ("E",):
            assert m[0] == len(m[1])

    # Keep everyone but one
    T = T_bak.copy()
    T = tx._get_sequence_budgets(T,3)
    merge_blocks = tx._identify_merge_blocks(T)
    expected_merges = {("A","B",):1,("C","D","H",):1,("E","F",):1}
    for m in merge_blocks:
        leaves = [leaf.name for leaf in m[2].get_leaves()]
        leaves.sort()
        leaves = tuple(leaves)
        assert m[0] == expected_merges[leaves]

def test_taxonomic_sample(for_real_inference,tmpdir):

    df = topiary.read_dataframe(for_real_inference["small-pre-redundancy.csv"])

    # Runs with default
    out = tx.taxonomic_sample(df)

    # no easy way to check that this *works* but at least this should not die
    # when other strings will
    out = tx.taxonomic_sample(df,paralog_column="nickname")

    # bad paralog_column values
    bad_paralog_column = ["not_a_column",None,14.2,str]
    for b in bad_paralog_column:
        print(b)
        with pytest.raises(ValueError):
            tx.taxonomic_sample(df,paralog_column=b)

    out = tx.taxonomic_sample(df,target_seq_number=1)
    assert np.sum(out.keep) == 1

    out = tx.taxonomic_sample(df,target_seq_number=5)
    assert np.sum(out.keep) == 5

    out = tx.taxonomic_sample(df,target_seq_number=9)
    assert np.sum(out.keep) == 9

    bad_target_seq_number = [-1,0,"test",int,14.2]
    for b in bad_target_seq_number:
        print(b)
        with pytest.raises(ValueError):
            tx.taxonomic_sample(df,target_seq_number=b)

    out = tx.taxonomic_sample(df,even_paralog_split=True,target_seq_number=5)
    assert np.sum(out.keep) == 5

    out = tx.taxonomic_sample(df,even_paralog_split=False,target_seq_number=5)
    assert np.sum(out.keep) == 5

    bad_even_paralog_split = [[],"test",int,14.2]
    for b in bad_even_paralog_split:
        print("sending paralog split",b,flush=True)
        with pytest.raises(ValueError):
            tx.taxonomic_sample(df,even_paralog_split=b)

    # Get rid of everything within species
    out = tx.taxonomic_sample(df,within_species_redundancy_cutoff=0.1,target_seq_number=1000)
    assert np.sum(out.keep) == 11

    # Don't get rid of anything within species
    out = tx.taxonomic_sample(df,within_species_redundancy_cutoff=1.0,target_seq_number=1000)
    assert np.sum(out.keep) == 19

    bad_redund = [-1,2,"stupid",None]
    for b in bad_redund:
        print("sending bad redundancy",b)
        with pytest.raises(ValueError):
            tx.taxonomic_sample(df,within_species_redundancy_cutoff=b)

    out = tx.taxonomic_sample(df,verbose=True)
    out = tx.taxonomic_sample(df,verbose=False)

    bad_verbose = [[],"test",int,14.2]
    for b in bad_verbose:
        print("sending verbose",b,flush=True)
        with pytest.raises(ValueError):
            tx.taxonomic_sample(df,verbose=b)
