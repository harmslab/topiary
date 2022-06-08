
import topiary
from topiary._arg_processors import process_bool, process_float, process_int, process_iter

import pytest
import numpy as np
import pandas as pd

def test_process_bool():

    true_values = [True,1.0,1,10,-1]
    for t in true_values:
        assert process_bool(t)

    false_values = [False,0.0,0]
    for f in false_values:
        assert not process_bool(f)

    bad_value = [None,"stupid",[1.0,1.0],np.array([1.0,1.0]),{},float]
    for b in bad_value:
        with pytest.raises(ValueError):
            value = process_float(b)


def test_process_float():

    value = process_float(1.0)
    assert value == 1.0

    bad_value = [None,"stupid",[1.0,1.0],np.array([1.0,1.0]),{},float,np.nan]
    for b in bad_value:
        with pytest.raises(ValueError):
            value = process_float(b)

    good_value = [-np.inf,np.inf,0,1,"1.0",np.array([1.0])]
    for g in bad_value:
        with pytest.raises(ValueError):
            value = process_float(g)

    with pytest.raises(ValueError):
        process_float(1.0,minimum_allowed=2.0)

    with pytest.raises(ValueError):
        process_float(1.0,minimum_allowed=1.0,minimum_inclusive=False)

    value = process_float(1.0,minimum_allowed=1.0,minimum_inclusive=True)
    assert value == 1

    with pytest.raises(ValueError):
        process_float(1.0,maximum_allowed=0.5)

    with pytest.raises(ValueError):
        process_float(1.0,maximum_allowed=1.0,maximum_inclusive=False)

    value = process_float(1.0,minimum_allowed=1.0,maximum_inclusive=True)
    assert value == 1


def test_process_int():

    value = process_int(1)
    assert value == 1

    bad_value = [None,"stupid",[1.0,1.0],np.array([1.0,1.0]),{},float,int,np.inf,np.nan,1.3]
    for b in bad_value:
        print(b)
        with pytest.raises(ValueError):
            value = process_int(b)

    good_value = [-10,0,10,"10",10.0,"10.0"]
    for g in bad_value:
        with pytest.raises(ValueError):
            value = process_int(g)

    with pytest.raises(ValueError):
        process_int(1,minimum_allowed=2.0)

    with pytest.raises(ValueError):
        process_int(1,minimum_allowed=1,minimum_inclusive=False)

    value = process_int(1,minimum_allowed=1,minimum_inclusive=True)
    assert value == 1

    with pytest.raises(ValueError):
        process_int(1,maximum_allowed=0)

    with pytest.raises(ValueError):
        process_int(1,maximum_allowed=1,maximum_inclusive=False)

    value = process_int(1,minimum_allowed=1,maximum_inclusive=True)
    assert value == 1

def test_process_iter():

    value = process_iter([1])
    assert np.array_equal(value,[1])

    bad_value = [None,0,list,float,int,np.inf,np.nan,1.3]
    for b in bad_value:
        print(b)
        with pytest.raises(ValueError):
            value = process_iter(b)


    df = pd.DataFrame({"test":[1,2,3]})
    good_value = [[],(),"test",{},df,np.arange(10)]
    for g in bad_value:
        with pytest.raises(ValueError):
            value = process_iter(g)

    value = process_iter([1,2,3],required_iter_type=list)
    for t in [str,tuple,type(df),type(np.arange(1))]:
        with pytest.raises(ValueError):
            value = process_iter([1,2,3],required_iter_type=t)

    value = process_iter(tuple([1,2,3]),required_iter_type=tuple)
    for t in [str,list,type(df),type(np.arange(1))]:
        with pytest.raises(ValueError):
            value = process_iter(tuple([1,2,3]),required_iter_type=t)

    value = process_iter(np.arange(5),required_iter_type=type(np.arange(1)))
    for t in [str,tuple,type(df),list]:
        with pytest.raises(ValueError):
            value = process_iter(np.arange(5),required_iter_type=t)

    value = process_iter(df,required_iter_type=type(df))
    for t in [str,tuple,type(np.arange(1)),list]:
        with pytest.raises(ValueError):
            value = process_iter(df,required_iter_type=t)

    # Check required value type check
    value = process_iter([1,2,3],required_value_type=int)
    for v in [["test"],[1.0],[None],[1,1.0]]:
        with pytest.raises(ValueError):
            value = process_iter(v,required_value_type=int)

    # Check limits checks
    with pytest.raises(ValueError):
        process_iter([1],minimum_allowed=2)

    with pytest.raises(ValueError):
        process_iter([1],minimum_allowed=1,minimum_inclusive=False)

    value = process_iter([1],minimum_allowed=1,minimum_inclusive=True)

    with pytest.raises(ValueError):
        process_iter([1,2],maximum_allowed=1)

    with pytest.raises(ValueError):
        process_iter([1,2],maximum_allowed=2,maximum_inclusive=False)

    value = process_iter([1,2],maximum_allowed=2,maximum_inclusive=True)

    # check is_not_type
    process_iter([1],is_not_type=str)
    with pytest.raises(ValueError):
        process_iter("test",is_not_type=str)

    with pytest.raises(ValueError):
        process_iter("test",is_not_type=[str,list])

    with pytest.raises(ValueError):
        process_iter([1,2],is_not_type=[str,list])

    process_iter(tuple([1,2]),is_not_type=[str,list])

    # with pytest.raises(ValueError):
    #     process_iter(1,maximum_allowed=0)
    #
    # with pytest.raises(ValueError):
    #     process_iter(1,maximum_allowed=1,maximum_inclusive=False)
    #
    # value = process_iter(1,minimum_allowed=1,maximum_inclusive=True)
    # assert value == 1
        #
        #
        #
        # value,
        #  variable_name=None,
        #  required_iter_type=None,
        #  required_value_type=None,
        #  minimum_length=None,
        #  maximum_length=None,
        #  minimum_inclusive=True,
        #  maximum_inclusive=True,
        #  is_not_type=None):
