
import pytest

import topiary
from topiary.check import check_paralog_patterns

import numpy as np
import pandas as pd
import re

def test_check_paralog_patterns():

    # Send in bad paralog_patterns data types (should be dict)
    bad_inputs = [1,-1,1.5,False,pd.DataFrame]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            patterns = check_paralog_patterns(paralog_patterns=b)

    # Send bad paralog_patterns keys (should be bad)
    bad_inputs = [(1,2),-1,False]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            patterns = check_paralog_patterns(paralog_patterns=b)

    # Send bad paralog_patterns values
    bad_inputs = [(1,2),-1,False]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            patterns = check_paralog_patterns(paralog_patterns=b)

    # Various paralog_patterns calls hould work
    patterns = check_paralog_patterns(paralog_patterns={"test":"this"})
    assert len(patterns) == 1
    assert patterns[0][0].search("string matches this") is not None
    assert patterns[0][0].search("string does not match") is None
    assert patterns[0][1] == "test"

    patterns = check_paralog_patterns(paralog_patterns={"test":["this","other"]})
    assert len(patterns) == 1
    assert patterns[0][0].search("string matches this") is not None
    assert patterns[0][0].search("string matches other") is not None
    assert patterns[0][0].search("string does not match") is None
    assert patterns[0][1] == "test"

    patterns = check_paralog_patterns(paralog_patterns={"test":("this","other")})
    assert len(patterns) == 1
    assert patterns[0][0].search("string matches this") is not None
    assert patterns[0][0].search("string matches other") is not None
    assert patterns[0][0].search("string does not match") is None
    assert patterns[0][1] == "test"

    patterns = check_paralog_patterns(paralog_patterns={"test":(re.compile("this"),"other")})
    assert len(patterns) == 1
    assert patterns[0][0].search("string matches this") is not None
    assert patterns[0][0].search("string matches other") is not None
    assert patterns[0][0].search("string does not match") is None
    assert patterns[0][1] == "test"

    patterns = check_paralog_patterns(paralog_patterns={"test":"this",
                                                          "thing":"other"})
    assert len(patterns) == 2
    assert patterns[0][0].search("string matches this") is not None
    assert patterns[0][0].search("string matches other") is None
    assert patterns[0][1] == "test"
    assert patterns[1][0].search("string matches other") is not None
    assert patterns[1][0].search("string does not match") is None
    assert patterns[1][1] == "thing"

    patterns = check_paralog_patterns(paralog_patterns={"test":"this"})
    assert patterns[0][0].flags == re.compile("a",flags=re.IGNORECASE).flags

    patterns = check_paralog_patterns(paralog_patterns={"test":"this"},
                                              ignorecase=True)
    assert patterns[0][0].flags == re.compile("a",flags=re.IGNORECASE).flags

    patterns = check_paralog_patterns(paralog_patterns={"test":"this"},
                                              ignorecase=False)
    assert patterns[0][0].flags == re.compile("a").flags
