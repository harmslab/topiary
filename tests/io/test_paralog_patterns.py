import pytest

from topiary.io.paralog_patterns import _get_alias_regex
from topiary.io.paralog_patterns import _build_alias_regex
from topiary.io.paralog_patterns import _preprocess_alias_input
from topiary.io.paralog_patterns import load_paralog_patterns

import numpy as np
import pandas as pd
import re

import copy


def test__get_alias_regex():

    test_strings = {"ABC": "abc",
                    "ABC1" :"abc[\ \-_\.]*1",
                    "1ABc" :"1[\ \-_\.]*abc",
                    " ABC ":"abc",
                    "AB C" :"ab[\ \-_\.]*c",
                    "(ABC)":"\(abc\)",
                    "[AB C]" :"\[ab[\ \-_\.]*c\]"}
    spacers = [" ","-","_","."]
    for t in test_strings:
        assert _get_alias_regex(t,spacers) == test_strings[t]


    test_strings = {"AbC": "abc",
                    "ABC1" :"abc[,]*1",
                    "1ABC" :"1[,]*abc",
                    " aBC ":"abc",
                    "AB C" :"ab\\ c"}
    spacers = [","]
    for t in test_strings:
        assert _get_alias_regex(t,spacers) == test_strings[t]

    test_strings = {"A" :"a",
                    "1" :"1",
                    "1A":"1[\ \-_\.]*a",
                    "A1":"a[\ \-_\.]*1",
                    " A1 " :"a[\ \-_\.]*1"}
    spacers = [" ","-","_","."]
    for t in test_strings:
        assert _get_alias_regex(t,spacers) == test_strings[t]


def test__build_alias_regex():

    alias_dict = {"test":["A"],
                  "this":["B"]}

    paralog_patterns = _build_alias_regex(alias_dict)
    assert isinstance(paralog_patterns,dict)
    assert len(paralog_patterns) == 2
    assert paralog_patterns["test"].pattern == "a|test"
    assert paralog_patterns["this"].pattern == "b|this"

    # Same alias in both
    alias_dict = {"test":["A"],
                  "this":["A"]}
    with pytest.raises(ValueError):
        paralog_patterns = _build_alias_regex(alias_dict)


    # Test alias that requires a negative regex to resolve.
    alias_dict = {"test":["A"],
                  "this":["AX"]}
    paralog_patterns = _build_alias_regex(alias_dict)
    assert paralog_patterns["test"].pattern == '^(?!.*(ax)).*(a|test)'
    assert paralog_patterns["this"].pattern == 'ax|this'

    # Make sure order doesn't matter for negative
    alias_dict = {"test":["AX"],
                  "this":["A"]}
    paralog_patterns = _build_alias_regex(alias_dict)
    assert paralog_patterns["test"].pattern == "ax|test"
    assert paralog_patterns["this"].pattern == '^(?!.*(ax)).*(a|this)'

    # Multiple inputs
    alias_dict = {"test":["A","B","C"],
                  "this":["D"]}
    paralog_patterns = _build_alias_regex(alias_dict)
    assert paralog_patterns["test"].pattern == "a|b|c|test"
    assert paralog_patterns["this"].pattern == "d|this"

    # Make sure regex kwargs are going in
    patterns = _build_alias_regex(alias_dict={"test":["this"]})
    assert patterns["test"].flags == re.compile("a",flags=re.IGNORECASE).flags

    patterns = _build_alias_regex(alias_dict={"test":["this"]},
                                              ignorecase=True)
    assert patterns["test"].flags == re.compile("a",flags=re.IGNORECASE).flags

    patterns = _build_alias_regex(alias_dict={"test":["this"]},
                                              ignorecase=False)
    assert patterns["test"].flags == re.compile("a").flags

    # First should NOT ignorecase (should not have been recompiled); second
    # should ignore case (would have to be compiled)
    px = re.compile("a")
    patterns = _build_alias_regex(alias_dict={"X":px,
                                              "Y":["b"]},
                                  ignorecase=True)
    assert patterns["X"].flags == re.compile("a").flags
    assert patterns["Y"].flags == re.compile("a",flags=re.IGNORECASE).flags


def test__preprocess_alias_input():

    # Send in bad paralog_patterns data types (should be dict)
    bad_inputs = [1,-1,1.5,False,pd.DataFrame]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            patterns = _preprocess_alias_input(alias_dict=b)

    # Send bad paralog_patterns keys (should be bad)
    bad_inputs = [(1,2),-1,False]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            patterns = _preprocess_alias_input(alias_dict=b)

    # Send bad paralog_patterns values
    bad_inputs = [(1,2),-1,False]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            patterns = _preprocess_alias_input(alias_dict=b)

    # Various paralog_patterns calls should work
    patterns = _preprocess_alias_input(alias_dict={"test":"this"})
    assert len(patterns) == 1
    assert len(patterns["test"]) == 1
    assert patterns["test"][0] == "this"

    patterns = _preprocess_alias_input(alias_dict={"test":["this","other"]})
    assert len(patterns) == 1
    assert len(patterns["test"]) == 2
    assert patterns["test"][0] == "this"
    assert patterns["test"][1] == "other"

    patterns = _preprocess_alias_input(alias_dict={"test":("this","other")})
    assert len(patterns) == 1
    assert len(patterns["test"]) == 2
    assert patterns["test"][0] == "this"
    assert patterns["test"][1] == "other"

    px = re.compile("other")
    patterns = _preprocess_alias_input(alias_dict={"test":px})
    assert len(patterns) == 1
    assert patterns["test"] == px

    patterns = _preprocess_alias_input(alias_dict={"test":(re.compile("this"),"other")})
    assert len(patterns) == 1
    assert len(patterns["test"]) == 2
    assert patterns["test"][0] == "this"
    assert patterns["test"][1] == "other"

    patterns = _preprocess_alias_input(alias_dict={"test":"this",
                                                   "thing":"other"})
    assert len(patterns) == 2
    assert len(patterns["test"]) == 1
    assert len(patterns["thing"]) == 1
    assert patterns["test"][0] == "this"
    assert patterns["thing"][0] == "other"

def test_load_paralog_patterns():
    
    default_kwargs = {"spacers":[" ","-","_","."],
                      "ignorecase":True,
                      "re_flags":None}

    kwargs = copy.deepcopy(default_kwargs)


    # Send in same pattern four equivalent ways and make sure it
    # ends up the same
    out = load_paralog_patterns(alias_dict={"stupid":"1"},**kwargs)
    patterns1 = copy.deepcopy(out)
    out = load_paralog_patterns(alias_dict={"stupid":["1"]},**kwargs)
    patterns2 = copy.deepcopy(out)
    out = load_paralog_patterns(alias_dict={"stupid":re.compile("1|stupid")},**kwargs)
    patterns3 = copy.deepcopy(out)
    out = load_paralog_patterns(alias_dict={"stupid":[re.compile("1")]},**kwargs)
    patterns4 = copy.deepcopy(out)

    pattern_list = [patterns1,patterns2,patterns3,patterns4]
    for i in range(len(pattern_list)):
        assert len(pattern_list[i]) == 1
        for j in range(i+1,len(pattern_list)):
            assert pattern_list[i]["stupid"].pattern == pattern_list[j]["stupid"].pattern

    # Send in all sorts of combinations of good paralog_patterns and make sure
    # it doesn't choke
    good_paralog_patterns = [{"stupid":"a","is":"b"},
                             {"stupid":["a"],"is":"b"},
                             {"stupid":["a"],"is":["b"]},
                             {"stupid":("a",),"is":("b","c")},
                             {"stupid":re.compile("a"),"is":"b"},
                             {"stupid":re.compile("a"),"is":re.compile("b")},
                             {"stupid":[re.compile("a")],"is":[re.compile("b")]}]
    for g in good_paralog_patterns:
        out = load_paralog_patterns(alias_dict=g,**kwargs)

    # Stupid things to pass in
    bad_paralog_patterns = [{"stupid":[]},
                            {"stupid":[],"is":[]},
                            {"stupid":"a","is":[]},
                            None,[],dict,pd.DataFrame({"test":[1,2,3]}),1,0,-1]
    for b in bad_paralog_patterns:
        print(f"passing {b} to _prepare_for_recip_blast")
        with pytest.raises(ValueError):
            out = load_paralog_patterns(alias_dict=b,**kwargs)

    # Slightly mangled things to pass in
    bad_paralog_patterns = [{1:["pattern"]},
                            {("1","2"):["pattern"]},
                            {"1":1}]
    for b in bad_paralog_patterns:
        with pytest.raises(ValueError):
            out = load_paralog_patterns(alias_dict=b,**kwargs)
