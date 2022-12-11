
import pytest

import topiary
from topiary._private import check
from topiary.io.seed import _get_alias_regex
from topiary.io.seed import _build_alias_regex
from topiary.io.seed import read_seed
from topiary.io.seed import df_from_seed

import numpy as np
import pandas as pd

import re

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


def test_read_seed(seed_dataframes,user_seed_dataframes):

    def _validate_output(out):

        df = out[0]
        key_species = out[1]
        paralog_patterns = out[2]
        species_aware = out[3]

        # check output dataframe
        assert len(df) == 8

        names = list(np.unique(df.name))
        names.sort()
        assert np.array_equal(names,("LY86","LY96"))

        species = list(np.unique(df.species))
        species.sort()
        assert np.array_equal(species,("Danio rerio",
                                       "Gallus gallus",
                                       "Homo sapiens",
                                       "Mus musculus"))

        assert np.array_equal(df.always_keep,np.ones(8,dtype=bool))

        # Make sure it's a good topiary dataframe
        check.check_topiary_dataframe(df)

        # Test key species
        assert np.array_equal(key_species,("Danio rerio",
                                           "Gallus gallus",
                                           "Homo sapiens",
                                           "Mus musculus"))

        # Test paralog patterns
        expected_patterns = {
            "LY96":"esop[\ \-_\.]*1|ly[\ \-_\.]*96|lymphocyte[\ \-_\.]*antigen[\ \-_\.]*96|myeloid[\ \-_\.]*differentiation[\ \-_\.]*protein[\ \-_\.]*2",
            "LY86":"ly[\ \-_\.]*86|lymphocyte[\ \-_\.]*antigen[\ \-_\.]*86|md[\ \-_\.]*1|mmd[\ \-_\.]*1|rp[\ \-_\.]*105[\ \-_\.]*associated[\ \-_\.]*3"
        }
        for k in paralog_patterns:
            print(paralog_patterns[k].pattern)
            assert paralog_patterns[k].pattern == expected_patterns[k]

        assert species_aware is True

    # csv
    df_file = seed_dataframes["good-seed-df.csv"]
    out = read_seed(df_file)
    _validate_output(out)

    # Make sure we can read the dataframe as a dataframe object
    df_as_df = pd.read_csv(df_file)
    out = read_seed(df_as_df)
    _validate_output(out)

    # tsv
    df_file = seed_dataframes["good-seed-df.tsv"]
    out = read_seed(df_file)
    _validate_output(out)

    # xlsx
    df_file = seed_dataframes["good-seed-df.xlsx"]
    out = read_seed(df_file)
    _validate_output(out)

    # csv with .txt extension
    df_file = seed_dataframes["good-seed-df.txt"]
    out = read_seed(df_file)
    _validate_output(out)

    # bad df passes
    with pytest.raises(FileNotFoundError):
        read_seed("not-really-a-file")

    bad_df = [pd.DataFrame,None,1,[],{},1.5,str,int,float]
    for b in bad_df:
        print(f"passing bad_df {b}")
        with pytest.raises(ValueError):
            read_seed(b)

    # required column check
    good_df = pd.read_csv(seed_dataframes["good-seed-df.csv"])
    bad_df = good_df.drop(columns=["species"])
    with pytest.raises(ValueError):
        read_seed(bad_df)

    # bad species
    bad_df = good_df.copy()
    bad_df.loc[:,"species"] = "not a species"
    with pytest.raises(ValueError):
        with pytest.warns():
            read_seed(bad_df)

    # species that is findable but not resolvable
    bad_df = good_df.copy()
    bad_df.loc[:,"species"] = "Bos indicus x Bos taurus"
    with pytest.raises(ValueError):
        read_seed(bad_df)

    # species that is findable but not resolvable
    bad_df = seed_dataframes["duplicated-alias.xlsx"]
    with pytest.raises(ValueError):
        read_seed(bad_df)

    # Read in a collection of user-generated dataframes, checking for expected
    # outputs
    
    is_species_aware = {'snase.xlsx':False,
                        'rnaseh.xlsx':False,
                        'chs.xlsx':True,
                        'iapp-cgrp.xlsx':True,
                        'gproteins.xlsx':False,
                        'ly86ly96.xlsx':True,
                        'zo1.xlsx': True,
                        's100a5-a6.xlsx':True}

    for s in user_seed_dataframes:

        print(f"Testing read of {s}")

        df, key_species, paralog_patterns, species_aware = topiary.io.read_seed(user_seed_dataframes[s])
        check_read = pd.read_excel(user_seed_dataframes[s])
        check_read.columns = [c.lower().strip() for c in check_read.columns]

        assert len(df) == len(check_read)

        for idx in df.index:
            assert df.loc[idx,"species"] == check_read.loc[idx,"species"].strip()
            assert df.loc[idx,"name"] == check_read.loc[idx,"name"].strip()

            this_seq = "".join(check_read.loc[idx,"sequence"].split())
            assert df.loc[idx,"sequence"] == this_seq

        assert np.array_equal(df.keep,np.ones(len(df),dtype=bool))
        assert np.array_equal(df.key_species,np.ones(len(df),dtype=bool))
        assert np.array_equal(df.always_keep,np.ones(len(df),dtype=bool))
        assert len(df.uid) == len(np.unique(df.uid))

        # Make sure the code is correctly identifying whether to treat this 
        # dataset as species aware
        assert species_aware is is_species_aware[s]



def test_df_from_seed():

    pass
