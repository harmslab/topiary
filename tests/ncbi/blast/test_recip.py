import pytest
from conftest import get_public_param_defaults

from topiary.ncbi.blast.recip import recip_blast
from topiary.ncbi.blast.recip import _calc_hit_post_prob
from topiary.ncbi.blast.recip import _prepare_for_blast
from topiary.ncbi.blast.recip import _run_blast
from topiary.ncbi.blast.recip import recip_blast
from topiary.ncbi.blast.recip import _make_recip_blast_calls

import numpy as np
import pandas as pd
import copy, re

def test__prepare_for_blast(test_dataframes):


    default_kwargs = get_public_param_defaults(recip_blast,
                                               _prepare_for_blast)

    df = test_dataframes["good-df"]
    paralog_patterns = {"LY96":["lymphocyte antigen 96","esop1"],
                        "LY86":re.compile("lymphocyte antigen 86")}

    kwargs = copy.deepcopy(default_kwargs)

    # Should fail, no blast db specified
    with pytest.raises(ValueError):
        out = _prepare_for_blast(df,paralog_patterns,**kwargs)

    # Should work (one blast db specified)
    kwargs["local_blast_db"] = "local"
    out = _prepare_for_blast(df,paralog_patterns,**kwargs)

    # Should fail (two blast db specified)
    kwargs["ncbi_blast_db"] = "nr"
    with pytest.raises(ValueError):
        out = _prepare_for_blast(df,paralog_patterns,**kwargs)

    # Should work (one blast db specified)
    kwargs["local_blast_db"] = None
    out = _prepare_for_blast(df,paralog_patterns,**kwargs)

    # Now hack default kwargs so it should run without modification
    default_kwargs["df"] = df.copy()
    default_kwargs["paralog_patterns"] = copy.deepcopy(paralog_patterns)
    default_kwargs["local_blast_db"] = "local"

    # Now start tests for each argument

    # ------------------------------------------------------------------------
    # df

    kwargs = copy.deepcopy(default_kwargs)

    # Should return a copy, but without modifying
    out = _prepare_for_blast(**kwargs)
    assert out[0] is not kwargs["df"]
    assert (out[0] == kwargs["df"].loc[:,out[0].columns]).all().all()

    non_topiary_df = pd.DataFrame({"test":[1,2,3]})

    bad_inputs = [1,-1,1.5,None,False,pd.DataFrame,non_topiary_df]
    for b in bad_inputs:
        kwargs["df"] = b
        with pytest.raises(ValueError):
            out = _prepare_for_blast(**kwargs)

    # ------------------------------------------------------------------------
    # paralog_patterns

    kwargs = copy.deepcopy(default_kwargs)
    out = _prepare_for_blast(**kwargs)

    # Make sure it's compiling regular expressions as expected
    patterns = out[2]
    expected_re = [("|".join([re.escape(x) for x in paralog_patterns["LY96"]]),"LY96"),
                   (paralog_patterns["LY86"].pattern,"LY86")]
    for i, p in enumerate(patterns):
        assert patterns[p].pattern == expected_re[i][0]
        assert p == expected_re[i][1]

    # Dump paralog_patterns so we can jam in different ones
    kwargs.pop("paralog_patterns")

    # Send in same pattern four equivalent ways and make sure it
    # ends up the same
    out = _prepare_for_blast(paralog_patterns={"stupid":"1"},**kwargs)
    patterns1 = out[2]
    out = _prepare_for_blast(paralog_patterns={"stupid":["1"]},**kwargs)
    patterns2 = out[2]
    out = _prepare_for_blast(paralog_patterns={"stupid":re.compile("1")},**kwargs)
    patterns3 = out[2]
    out = _prepare_for_blast(paralog_patterns={"stupid":[re.compile("1")]},**kwargs)
    patterns4 = out[2]

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
        out = _prepare_for_blast(paralog_patterns=g,**kwargs)

    # Stupid things to pass in
    bad_paralog_patterns = [{},
                      {"stupid":[]},
                      {"stupid":[],"is":[]},
                      {"stupid":"a","is":[]},
                      None,[],dict,pd.DataFrame({"test":[1,2,3]}),1,0,-1]
    for b in bad_paralog_patterns:
        print(f"passing {b} to _prepare_for_recip_blast")
        with pytest.raises(ValueError):
            out = _prepare_for_blast(paralog_patterns=b,**kwargs)

    # Slightly mangled things to pass in
    bad_paralog_patterns = [{1:["pattern"]},
                            {("1","2"):["pattern"]},
                            {"1":1}]
    for b in bad_paralog_patterns:
        with pytest.raises(ValueError):
            out = _prepare_for_blast(paralog_patterns=b,**kwargs)



    # ------------------------------------------------------------------------
    # local_blast_db and ncbi_blast_db already tested above

    # ------------------------------------------------------------------------
    # ignorecase
    kwargs = copy.deepcopy(default_kwargs)
    kwargs.pop("ignorecase")

    # Make sure True passes work
    good_trues = [True,1,np.ones(1,dtype=bool)[0]]
    for g in good_trues:
        print("passing good",g,type(g))
        out = _prepare_for_blast(ignorecase=g,**kwargs)

    # Make sure False passes work
    good_falses = [False,0,np.zeros(1,dtype=bool)[0]]
    for g in good_falses:
        out = _prepare_for_blast(ignorecase=g,**kwargs)

    bad_ignorecase = ["True","False",np.nan,{},[],None,type(True)]
    for b in bad_ignorecase:
        print(f"testing bad ignorecase value: {b}")
        with pytest.raises(ValueError):
            out = _prepare_for_blast(ignorecase=b,**kwargs)


    # ------------------------------------------------------------------------
    # min_call_prob
    kwargs = copy.deepcopy(default_kwargs)
    kwargs.pop("min_call_prob")
    good_args = [0.00000001,0.5,0.999999999]
    for g in good_args:
        out = _prepare_for_blast(min_call_prob=g,**kwargs)
        assert out[3] == g

    bad_args = [0,1,-1,1.1,None,"a",[],type(float),{},np.nan]
    for b in bad_args:
        print(f"testing bad min_call_prob value: {b}")
        with pytest.raises(ValueError):
            out = _prepare_for_blast(min_call_prob=b,**kwargs)

    # ------------------------------------------------------------------------
    # partition_temp
    kwargs = copy.deepcopy(default_kwargs)
    kwargs.pop("partition_temp")
    good_args = [0.00000001,0.5,1e10]
    for g in good_args:
        out = _prepare_for_blast(partition_temp=g,**kwargs)
        assert out[4] == g

    bad_args = [0,-1,None,"a",[],type(float),{},np.nan]
    for b in bad_args:
        print(f"testing bad partition_temp value: {b}")
        with pytest.raises(ValueError):
            out = _prepare_for_blast(partition_temp=b,**kwargs)

    # ------------------------------------------------------------------------
    # drop_combo_fx
    kwargs = copy.deepcopy(default_kwargs)
    kwargs.pop("drop_combo_fx")
    good_args = [0.00000001,0.5,1]
    for g in good_args:
        out = _prepare_for_blast(drop_combo_fx=g,**kwargs)
        assert out[5] == g

    bad_args = [-1,10000000000.0,None,"a",[],type(float),{},np.nan]
    for b in bad_args:
        print(f"testing bad drop_combo_fx value: {b}")
        with pytest.raises(ValueError):
            out = _prepare_for_blast(drop_combo_fx=b,**kwargs)

    # ------------------------------------------------------------------------
    # use_start_end
    kwargs = copy.deepcopy(default_kwargs)
    kwargs.pop("use_start_end")

    hack_df = kwargs["df"].copy()
    hack_df["start"] = 10
    kwargs["df"] = hack_df

    # Make sure True passes work
    good_trues = [True,1,np.ones(1,dtype=bool)[0]]
    for g in good_trues:
        out = _prepare_for_blast(use_start_end=g,**kwargs)
        assert len(out[1][0]) == 150 # should use start/stop, truncated seq

    # Make sure False passes work
    good_falses = [False,0,np.zeros(1,dtype=bool)[0]]
    for g in good_falses:
        out = _prepare_for_blast(use_start_end=g,**kwargs)
        assert len(out[1][0]) == 160 # should not use start/stop, full length seq

    bad_use_start_end = ["True","False",np.nan,{},[],None,type(True)]
    for b in bad_use_start_end:
        print(f"testing bad use_start_end value: {b}")
        with pytest.raises(ValueError):
            out = _prepare_for_blast(use_start_end=b,**kwargs)

    # ------------------------------------------------------------------------
    # deal with sequence_list

    kwargs = copy.deepcopy(default_kwargs)
    hack_df = kwargs["df"].copy()
    hack_df["start"] = [  0, 10, 20, 40, 50]
    hack_df["end"] =   [160,150,140,120,110]
    kwargs["df"] = hack_df
    kwargs.pop("use_start_end")

    out = _prepare_for_blast(use_start_end=True,**kwargs)
    sequence_list = out[1]
    assert len(sequence_list) == 5
    assert set([type(s) for s in sequence_list]) == set([str])
    assert np.array_equal(np.array([len(s) for s in sequence_list]),
                          np.array([160,150-10,140-20,120-40,110-50]))

    out = _prepare_for_blast(use_start_end=False,**kwargs)
    sequence_list = out[1]
    assert np.array_equal(np.array([len(s) for s in sequence_list]),
                          np.array([160,160,160,160,160]))

    # Drop start and stop, but ask to use_start_end. Should work, but give
    # sequences all same length
    kwargs["df"] = kwargs["df"].drop(columns=["start","end"])
    out = _prepare_for_blast(use_start_end=True,**kwargs)
    sequence_list = out[1]
    assert np.array_equal(np.array([len(s) for s in sequence_list]),
                          np.array([160,160,160,160,160]))

    # Make a terrible start/end combo
    kwargs = copy.deepcopy(default_kwargs)
    hack_df = kwargs["df"].copy()
    hack_df["start"] = [160, 10, 20, 40, 50]
    hack_df["end"] =   [-1,150,140,120,110]
    kwargs["df"] = hack_df
    kwargs.pop("use_start_end")
    with pytest.raises(ValueError):
        out = _prepare_for_blast(use_start_end=True,**kwargs)

    # Check for appropriate keep behavior. Should drop one sequence for recip
    # blast because keep == False
    kwargs = copy.deepcopy(default_kwargs)
    hack_df = kwargs["df"].copy()
    hack_df.loc[0,"keep"] = False
    kwargs["df"] = hack_df
    out = _prepare_for_blast(**kwargs)
    assert len(out[1])  == 4

def test__calc_hit_post_prob():

    patterns = {"A":re.compile("A",flags=re.IGNORECASE),
                "B":re.compile("B",flags=re.IGNORECASE)}

    # Both hits identical
    hits = pd.DataFrame({"bits":[1,1],
                         "hit_def":["A","B"],
                         "e_value":[5,10]})
    paralogs, pp, masks = _calc_hit_post_prob(hits,patterns,partition_temp=1)
    assert np.array_equal(pp,[0.5,0.5])

    # strongly favor A
    hits = pd.DataFrame({"bits":[100,1],
                         "hit_def":["A","B"],
                         "e_value":[5,10]})
    paralogs, pp, masks = _calc_hit_post_prob(hits,patterns,partition_temp=1)
    assert np.isclose(pp[0],1)
    assert np.isclose(pp[1],0)

    # Strongly favor A, but partition temperature so high it doesn't matter
    paralogs, pp, masks = _calc_hit_post_prob(hits,patterns,partition_temp=1e9)
    assert np.isclose(pp[0],0.5)
    assert np.isclose(pp[1],0.5)

    # only A hits - even with high temp favor A
    hits = pd.DataFrame({"bits":[1,1],
                         "hit_def":["A","A"],
                         "e_value":[5,10]})
    paralogs, pp, masks = _calc_hit_post_prob(hits,patterns,partition_temp=1e9)
    assert np.isclose(pp[0],1)
    assert np.isclose(pp[1],0)

    # No hits at all
    hits = pd.DataFrame({"bits":[1,1],
                         "hit_def":["X","Y"],
                         "e_value":[5,10]})
    paralogs, pp, masks = _calc_hit_post_prob(hits,patterns,partition_temp=1e9)
    assert np.isclose(pp[0],0)
    assert np.isclose(pp[1],0)


def test__make_recip_blast_calls(test_dataframes,recip_blast_hit_dfs):

    # --------------------------------------------------------------------------
    # Run _prepare_for_blast so we have all of the inputs we need for
    # _make_recip_blast_calls
    prep_kwargs = get_public_param_defaults(recip_blast,
                                            _prepare_for_blast)

    df = test_dataframes["good-df"]
    paralog_patterns = {"LY96":["lymphocyte antigen 96","esop1"],
                        "LY86":re.compile("lymphocyte antigen 86")}

    # Should work (one blast db specified, doesn't actually change output)
    prep_kwargs["local_blast_db"] = "local"
    prep_out = _prepare_for_blast(df,paralog_patterns,**prep_kwargs)

    # Okay, output from _prepare_for_blast
    df, sequence_list, paralog_patterns, min_call_prob, partition_temp, drop_combo_fx = prep_out

    # Get default arguments for _make_recip_blast_calls
    default_kwargs = get_public_param_defaults(recip_blast,
                                               _make_recip_blast_calls)
    default_kwargs["df"] = df.copy()
    default_kwargs["paralog_patterns"] = copy.deepcopy(paralog_patterns)

    # Go through the two blast types
    for blast_type in recip_blast_hit_dfs:

        # Does it run?
        hit_dfs = recip_blast_hit_dfs[blast_type]
        kwargs = copy.deepcopy(default_kwargs)
        if blast_type.startswith("ncbi"):
            kwargs["ncbi_blast_db"] = "nr"
        out_df = _make_recip_blast_calls(hit_dfs=hit_dfs,**kwargs)

        assert np.all(out_df["keep"])
        assert np.all(out_df["recip_found_paralog"])


    # Force a match to not occur
    for blast_type in recip_blast_hit_dfs:

        hit_dfs = recip_blast_hit_dfs[blast_type]
        kwargs = copy.deepcopy(default_kwargs)
        if blast_type.startswith("ncbi"):
            kwargs["ncbi_blast_db"] = "nr"

        # Set up first sequence so hit definition does not match
        hit_dfs[0]["hit_def"] = "NOT MATCHING"
        out_df = _make_recip_blast_calls(hit_dfs=hit_dfs,**kwargs)

        assert np.array_equal(np.array(out_df["keep"]),
                              np.array([False,True,True,True,True]))
        assert np.array_equal(np.array(out_df["recip_found_paralog"]),
                              np.array([False,True,True,True,True]))

    # Force a match to not occur
    for blast_type in recip_blast_hit_dfs:

        hit_dfs = recip_blast_hit_dfs[blast_type]
        kwargs = copy.deepcopy(default_kwargs)
        if blast_type.startswith("ncbi"):
            kwargs["ncbi_blast_db"] = "nr"

        # Set up first sequence so hit definition does not match
        hit_dfs[0]["hit_def"] = "lymphocyte antigen 86"
        out_df = _make_recip_blast_calls(hit_dfs=hit_dfs,**kwargs)

        assert np.all(out_df["keep"])
        assert np.all(out_df["recip_found_paralog"])
        assert np.array_equal(np.array(out_df["recip_paralog"]),
                              np.array(["LY86","LY96","LY96","LY96","LY96"]))

def test__run_blast():

    pass

def test_recip_blast():

    pass
