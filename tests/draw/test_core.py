
import pytest

import topiary
from topiary.draw.core import ete3_to_toytree
from topiary.draw.core import construct_colormap, construct_sizemap
from topiary.draw.core import create_name_dict
from topiary.draw.core import _protect_name, _deprotect_name, color_to_css
from topiary.draw.core import get_round_to
from topiary.draw.core import load_trees
from topiary.draw.core import map_tree_to_tree
from topiary.draw.core import parse_position_string
from topiary.draw.core import parse_span_color

import ete3
import toyplot

import numpy as np
import pandas as pd
import png  # toyplot uses png, so this is not a new dependency

import string, os

def test__protect_name():

    assert _protect_name("test") == "'test'"
    assert _protect_name("te st") == "'te%20st'"
    assert _protect_name(" te st ") == "'%20te%20st%20'"
    assert _protect_name("te,st") == "'te,st'"

def test__deprotect_name():

    assert _deprotect_name("'test'") == "test"
    assert _deprotect_name("'te%20st'") == "te st"
    assert _deprotect_name("'%20te%20st%20'") == " te st "
    assert _deprotect_name("'te,st'") == "te,st"

def test_color_to_css():

    assert color_to_css('rgba(100.0%,0.0%,0.0%,1.000)') == 'rgba(100.0%,0.0%,0.0%,1.000)'
    assert color_to_css([1,0,0]) == 'rgba(100.0%,0.0%,0.0%,1.000)'
    assert color_to_css([1,0,0,1]) == 'rgba(100.0%,0.0%,0.0%,1.000)'
    assert color_to_css("red") == 'rgba(100.0%,0.0%,0.0%,1.000)'

    with pytest.raises(ValueError):
        color_to_css([1,0,0,1,1,1])

def test_get_round_to():

    assert get_round_to(1e50,total_requested=3) == 0
    assert get_round_to(0.1,total_requested=3) == 1
    assert get_round_to(0.01,total_requested=3) == 2
    assert get_round_to(0.001,total_requested=3) == 3
    assert get_round_to(0.0011,total_requested=3) == 3
    assert get_round_to(0.0001,total_requested=3) == 4
    assert get_round_to(0.00001,total_requested=3) == 5
    assert get_round_to(1e-5,total_requested=3) == 5

    assert get_round_to(1.1,total_requested=3) == 1
    assert get_round_to(1.10,total_requested=3) == 1
    assert get_round_to(1.12,total_requested=3) == 2
    assert get_round_to(1.12,total_requested=2) == 1
    assert get_round_to(1.12,total_requested=1) == 0

    assert get_round_to(-1e50,total_requested=3) == 0
    assert get_round_to(-0.1,total_requested=3) == 1
    assert get_round_to(-0.01,total_requested=3) == 2
    assert get_round_to(-0.001,total_requested=3) == 3
    assert get_round_to(-0.0011,total_requested=3) == 3
    assert get_round_to(-0.0001,total_requested=3) == 4
    assert get_round_to(-0.00001,total_requested=3) == 5
    assert get_round_to(-1e-5,total_requested=3) == 5

    assert get_round_to(-1.1,total_requested=3) == 1
    assert get_round_to(-1.10,total_requested=3) == 1
    assert get_round_to(-1.12,total_requested=3) == 2
    assert get_round_to(-1.12,total_requested=2) == 1
    assert get_round_to(-1.12,total_requested=1) == 0

def test_ete3_to_toytree():

    # Should work; no interesting data loaded
    tree = "(((A,B),((C,D),E)),(F,G));"
    T = ete3.Tree(tree)
    tT = ete3_to_toytree(T)

    # Duplicate taxon name
    tree = "(((A,A),((C,D),E)),(F,G));"
    with pytest.raises(ValueError):
        T = ete3.Tree(tree)
        tT = ete3_to_toytree(T)

    # Make sure we're loading dist correctly
    tree = "(((A:1.0,B:4.0):1.0,((C:1.0,D:0.5):1.0,E:1.0):1.0):1.0,(F:1.0,G:1.0));"
    T = ete3.Tree(tree)
    tT = ete3_to_toytree(T)

    expected_dist = dict([(a,1.0) for a in string.ascii_uppercase[:7]])
    expected_dist["B"] = 4.0
    expected_dist["D"] = 0.5
    for n in tT.idx_dict:
        node = tT.idx_dict[n]
        if node.is_leaf():
            name = node.name
            assert expected_dist[name] == node.dist

    # Make sure we're loading dist and support correctly
    tree = "(((A:1.0,B:4.0)99:1.0,((C:1.0,D:0.5)100:1.0,E:1.0)75:1.0)80:1.0,(F:1.0,G:1.0)90)75;"
    T = ete3.Tree(tree)
    tT = ete3_to_toytree(T)

    expected_dist = dict([(a,1.0) for a in string.ascii_uppercase[:7]])
    expected_dist["B"] = 4.0
    expected_dist["D"] = 0.5
    expected_support = {"AB":99,"CD":100,"CDE":75,"ABCDE":80,"FG":90,"ABCDEFG":75}

    for n in tT.idx_dict:
        node = tT.idx_dict[n]
        if node.is_leaf():
            name = node.name
            assert expected_dist[name] == node.dist
        else:
            leaves = node.get_leaf_names()
            leaves.sort()
            key = "".join(leaves)
            assert node.support == expected_support[key]

    # Make sure we're loading dist and custom names correctly
    tree = "(((A:1.0,B:4.0)AB:1.0,((C:1.0,D:1.0)CD:1.0,E:1.0)CDE:1.0)ABCDE:1.0,(F:1.0,G:1.0)FG)ABCDEFG;"
    T = ete3.Tree(tree,format=1)
    tT = ete3_to_toytree(T)

    expected_dist = dict([(a,1.0) for a in string.ascii_uppercase[:7]])
    expected_dist["B"] = 4.0
    expected_dist["D"] = 0.5

    for n in tT.idx_dict:
        node = tT.idx_dict[n]
        if not node.is_leaf():
            leaves = node.get_leaf_names()
            leaves.sort()
            expected_name = "".join(leaves)
            assert node.name == expected_name

    # Make sure we're loading custom fields correctly
    tree = "(((A:1.0,B:4.0)99:1.0,((C:1.0,D:0.5)100:1.0,E:1.0)75:1.0)80:1.0,(F:1.0,G:1.0)90)75;"
    T = ete3.Tree(tree)
    for i, node in enumerate(T.traverse()):
        leaves = node.get_leaf_names()
        leaves.sort()
        node.add_feature("test_feature","".join(leaves))

    expected_dist = dict([(a,1.0) for a in string.ascii_uppercase[:7]])
    expected_dist["B"] = 4.0
    expected_dist["D"] = 0.5

    tT = ete3_to_toytree(T)
    for n in tT.idx_dict:

        node = tT.idx_dict[n]

        leaves = node.get_leaf_names()
        leaves.sort()
        key = "".join(leaves)

        if node.is_leaf():
            name = node.name
            assert expected_dist[name] == node.dist
        else:
            assert node.support == expected_support[key]

        assert node.test_feature == key

    # Make sure we can handel spaces in taxon names correctly
    tree = "((('A name':1.0,'B name':4.0):1.0,((C:1.0,D:0.5):1.0,E:1.0):1.0):1.0,(F:1.0,G:1.0));"
    T = ete3.Tree(tree)
    tT = ete3_to_toytree(T)


def test_construct_colormap():

    # only pass in value_max
    cmap, cmap_span =construct_colormap(color=["white","red"],prop=[0,100])
    assert cmap(100) == 'rgba(100.0%,0.0%,0.0%,1.000)'
    assert cmap(50) == 'rgba(100.0%,50.0%,50.0%,1.000)'
    assert cmap(0) == 'rgba(100.0%,100.0%,100.0%,1.000)'
    assert np.array_equal(cmap_span, [0,100])

    # value_max with custom colors
    cmap, cmap_span =construct_colormap(prop=(0,1),
                              color=("#FFFFFF","#FF0000"))

    assert cmap(1) == 'rgba(100.0%,0.0%,0.0%,1.000)'
    assert cmap(0.5) == 'rgba(100.0%,50.0%,50.0%,1.000)'
    assert cmap(0) == 'rgba(100.0%,100.0%,100.0%,1.000)'
    assert np.array_equal(cmap_span, [0,1])

    # value_max with flipped custom colors
    cmap, cmap_span =construct_colormap(prop=(0,1),
                              color=("#00FF00","#FFFFFF"))

    assert cmap(0) == 'rgba(0.0%,100.0%,0.0%,1.000)'
    assert cmap(0.5) == 'rgba(50.0%,100.0%,50.0%,1.000)'
    assert cmap(1) == 'rgba(100.0%,100.0%,100.0%,1.000)'
    assert np.array_equal(cmap_span, [0,1])

    # value_max and value_min with custom colors
    cmap, cmap_span =construct_colormap(prop=(-1,1),
                              color=("#FFFFFF","#00FF00"))

    assert cmap(1) == 'rgba(0.0%,100.0%,0.0%,1.000)'
    assert cmap(0) == 'rgba(50.0%,100.0%,50.0%,1.000)'
    assert cmap(-1) == 'rgba(100.0%,100.0%,100.0%,1.000)'
    assert np.array_equal(cmap_span, [-1,1])

    # bad value_span
    bad_values = ["test",[],float]
    for b in bad_values:
        print("passing",b)
        with pytest.raises(ValueError):
            construct_colormap(prop=(b,b),color=["red","white"])

    # pass in palette
    palette = toyplot.color.Palette(colors=["#00FF00","#FFFFFF"])
    cmap, cmap_span =construct_colormap(color=None,prop=(-1,1),palette=palette)
    assert cmap(-1) == 'rgba(0.0%,100.0%,0.0%,1.000)'
    assert cmap(0) == 'rgba(50.0%,100.0%,50.0%,1.000)'
    assert cmap(1) == 'rgba(100.0%,100.0%,100.0%,1.000)'
    assert np.array_equal(cmap_span, [-1,1])

    # pass in bad palette
    palette = "stupid"
    with pytest.raises(TypeError):
        cmap, cmap_span =construct_colormap(color=None,
                                  prop=(-1,1),
                                  palette=palette)

    # Not really a gradient -- just color_min
    cmap, cmap_span =construct_colormap(prop=(1,1),
                              color=("#00FF00","#FFFFFF"))

    assert cmap(-1) == 'rgba(0.0%,100.0%,0.0%,1.000)'
    assert cmap(0) == 'rgba(0.0%,100.0%,0.0%,1.000)'
    assert cmap(1) == 'rgba(0.0%,100.0%,0.0%,1.000)'
    assert len(cmap_span) == 0

    # single color
    cmap, cmap_span =construct_colormap(prop=(0,1),color="#00FF00")
    assert cmap(0) == 'rgba(0.0%,100.0%,0.0%,1.000)'
    assert cmap(1) == 'rgba(0.0%,100.0%,0.0%,1.000)'
    assert len(cmap_span) == 0

    # rgb
    cmap, cmap_span =construct_colormap(prop=(0,1),color=(1,1,1))
    assert cmap(0) == 'rgba(100.0%,100.0%,100.0%,1.000)'
    assert cmap(1) == 'rgba(100.0%,100.0%,100.0%,1.000)'
    assert len(cmap_span) == 0

    # rgba
    cmap, cmap_span =construct_colormap(prop=(0,1),color=(1,1,1,1))
    assert cmap(0) == 'rgba(100.0%,100.0%,100.0%,1.000)'
    assert cmap(1) == 'rgba(100.0%,100.0%,100.0%,1.000)'
    assert len(cmap_span) == 0

    # color dictionary
    cmap, cmap_span = construct_colormap(prop=("A","B"),color={"A":"white","B":"red"})
    assert cmap("A") == 'rgba(100.0%,100.0%,100.0%,1.000)'
    assert cmap("B") == 'rgba(100.0%,0.0%,0.0%,1.000)'
    with pytest.raises(KeyError):
        cmap(1)
    assert len(set(cmap_span).difference(set(["A","B"]))) == 0

    with pytest.raises(ValueError):
        construct_colormap(prop=("A","B","C"),color={"A":"white","B":"red"})


def test_construct_sizemap():

    sm, sm_span =construct_sizemap(prop=(0,1),size=(0,1))
    assert sm(0) == 0
    assert sm(0.5) == 0.5
    assert sm(1) == 1
    assert np.array_equal(sm_span,[0,1])

    sm, sm_span =construct_sizemap(prop=(0,1),size=(0,10))
    assert sm(0) == 0
    assert sm(0.5) == 5
    assert sm(1) == 10
    assert np.array_equal(sm_span,[0,1])

    sm, sm_span =construct_sizemap(prop=(-1,1),size=(0,10))
    assert sm(-1) == 0
    assert sm(0) == 5
    assert sm(1) == 10
    assert np.array_equal(sm_span,[-1,1])

    sm, sm_span =construct_sizemap(prop=(-1,1),size=(10,0))
    assert sm(-1) == 10
    assert sm(0) == 5
    assert sm(1) == 0
    assert np.array_equal(sm_span,[-1,1])

    sm, sm_span =construct_sizemap(prop=(1,1),size=(10,0))
    assert sm(-1) == 10
    assert sm(0) == 10
    assert sm(1) == 10
    assert len(sm_span) == 0

    sm, sm_span =construct_sizemap(prop=(-1,1),size=(10,10))
    assert sm(-1) == 10
    assert sm(0) == 10
    assert sm(1) == 10
    assert len(sm_span) == 0

    # Single value for size
    sm, sm_span =construct_sizemap(prop=(-1,1),size=5)
    assert sm(-1) == 5
    assert sm(0) == 5
    assert sm(1) == 5
    assert len(sm_span) == 0

    # Bad size
    with pytest.raises(ValueError):
        sm, sm_span =construct_sizemap(prop=(-1,1),size=-1)

    # Another bad size
    with pytest.raises(ValueError):
        sm, sm_span =construct_sizemap(prop=(-1,1),size=(-1,0))

    # size dictionary
    sm, sm_span =construct_sizemap(prop=(-1,1),size={-1:10,1:5})
    assert sm(-1) == 10
    assert sm(1) == 5
    with pytest.raises(KeyError):
        assert sm(0) == 5
    assert len(set(sm_span).difference(set([-1,1]))) == 0

    sm, sm_span =construct_sizemap(prop=("A","B"),size={"A":10,"B":5})
    assert sm("A") == 10
    assert sm("B") == 5
    assert len(set(sm_span).difference(set(["A","B"]))) == 0

    # Bad size given non-float property
    with pytest.raises(ValueError):
        sm, sm_span =construct_sizemap(prop=("A","B"),size=(0,1))

    # This should work -- size is fixed
    sm, sm_span =construct_sizemap(prop=("A","B"),size=1)
    assert len(sm_span) == 0


def test_create_name_dict(for_real_inference):

    df = topiary.read_dataframe(for_real_inference["small-pre-redundancy.csv"])

    # Make sure unique values come in
    name_dict = create_name_dict(df)
    assert issubclass(type(name_dict),dict)
    values = list(name_dict.values())
    assert len(values) == len(set(values))

    # Make all sequence names not unique -- should be smart enough to recognize
    # that and make unique
    test_df = df.copy()
    test_df.species = "Homo sapiens"
    name_dict = create_name_dict(test_df)
    assert issubclass(type(name_dict),dict)
    values = list(name_dict.values())
    assert len(values) == len(set(values))

    name_dict = create_name_dict(df,tip_columns=["uid"])

    # Make sure it's checking for bad value correctly
    with pytest.raises(ValueError):
        name_dict = create_name_dict(df,tip_columns=["uid"],separator=",")

def test_load_trees():

    pass

def test_map_tree_to_tree():

    pass

def test_parse_position_string():

    pass

def test_parse_span_color():

    pass
