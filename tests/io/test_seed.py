
import pytest

import topiary
from topiary import check
from topiary.io.seed import _get_string_variants, _get_alias_regex
from topiary.io.seed import load_seed_dataframe

import numpy as np
import pandas as pd

import re

def test__get_string_variants():

    test_strings = {"ABC": "abc",
                    "ABC1" :"abc[\ \-_\.]*1",
                    "1ABc" :"1[\ \-_\.]*abc",
                    " ABC ":"abc",
                    "AB C" :"ab[\ \-_\.]*c"}
    spacers = [" ","-","_","."]
    for t in test_strings:
        assert _get_string_variants(t,spacers) == test_strings[t]


    test_strings = {"AbC": "abc",
                    "ABC1" :"abc[,]*1",
                    "1ABC" :"1[,]*abc",
                    " aBC ":"abc",
                    "AB C" :"ab c"}
    spacers = [","]
    for t in test_strings:
        assert _get_string_variants(t,spacers) == test_strings[t]

    test_strings = {"A" :"a",
                    "1" :"1",
                    "1A":"1[\ \-_\.]*a",
                    "A1":"a[\ \-_\.]*1",
                    " A1 " :"a[\ \-_\.]*1"}
    spacers = [" ","-","_","."]
    for t in test_strings:
        assert _get_string_variants(t,spacers) == test_strings[t]

def test__get_alias_regex():

    test_names = {("LY96","ESOP1"):"esop[\ \-_\.]*1|ly[\ \-_\.]*96",
                  ("LY96",):   "ly[\ \-_\.]*96",
                  ("ESOP1",):"esop[\ \-_\.]*1",
                  ("1A","1A","A1"," A "):"1[\ \-_\.]*a|a|a[\ \-_\.]*1"}
    for t in test_names:
        observed = _get_alias_regex(t)
        expected = re.compile(test_names[t],flags=re.IGNORECASE)
        assert observed.pattern == expected.pattern

def test_load_seed_dataframe(seed_dataframes):

    def _validate_output(df, phylo_context, key_species, paralog_patterns):

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

        # Test phylo context
        assert phylo_context == "Vertebrates"

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
            assert paralog_patterns[k].pattern == expected_patterns[k]



    # csv
    df_file = seed_dataframes["good-seed-df.csv"]
    df, phylo_context, key_species, paralog_patterns = load_seed_dataframe(df_file)
    _validate_output(df, phylo_context, key_species, paralog_patterns)

    # Make sure we can read the dataframe as a dataframe object
    df_as_df = pd.read_csv(df_file)
    df, phylo_context, key_species, paralog_patterns = load_seed_dataframe(df_as_df)
    _validate_output(df, phylo_context, key_species, paralog_patterns)

    # tsv
    df_file = seed_dataframes["good-seed-df.tsv"]
    df, phylo_context, key_species, paralog_patterns = load_seed_dataframe(df_file)
    _validate_output(df, phylo_context, key_species, paralog_patterns)

    # xlsx
    df_file = seed_dataframes["good-seed-df.xlsx"]
    df, phylo_context, key_species, paralog_patterns = load_seed_dataframe(df_file)
    _validate_output(df, phylo_context, key_species, paralog_patterns)

    # csv with .txt extension
    df_file = seed_dataframes["good-seed-df.txt"]
    df, phylo_context, key_species, paralog_patterns = load_seed_dataframe(df_file)
    _validate_output(df, phylo_context, key_species, paralog_patterns)

    # bad df passes
    with pytest.raises(FileNotFoundError):
        load_seed_dataframe("not-really-a-file")

    bad_df = [pd.DataFrame,None,1,[],{},1.5,str,int,float]
    for b in bad_df:
        print(f"passing bad_df {b}")
        with pytest.raises(ValueError):
            load_seed_dataframe(b)

    # required column check
    good_df = pd.read_csv(seed_dataframes["good-seed-df.csv"])
    bad_df = good_df.drop(columns=["species"])
    with pytest.raises(ValueError):
        load_seed_dataframe(bad_df)

    # bad species
    bad_df = good_df.copy()
    bad_df.loc[:,"species"] = "not a species"
    with pytest.raises(ValueError):
        load_seed_dataframe(bad_df)

    # Change up species being sent in and make sure phylo_context changes
    hacked_df = good_df.copy()
    hacked_df.loc[0,"species"] = "Escherichia coli"
    df, phylo_context, key_species, paralog_patterns = load_seed_dataframe(hacked_df)
    assert phylo_context == "All life"

    hacked_df.loc[:,"species"] = "Homo sapiens"
    df, phylo_context, key_species, paralog_patterns = load_seed_dataframe(hacked_df)
    assert phylo_context == "Mammals"
