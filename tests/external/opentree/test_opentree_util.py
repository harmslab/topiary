

import pytest
import topiary
from topiary import opentree

import pandas as pd
import numpy as np

import re


def test_get_phylo_context():

    species_list = ["Homo sapiens","Gallus gallus","Danio rerio","Thermus thermophilus"]
    result = opentree.util.get_phylo_context(species_list)
    assert result == "All life"

    species_list = ["Homo sapiens","Gallus gallus","Danio rerio"]
    result = opentree.util.get_phylo_context(species_list)
    assert result == "Vertebrates"

    species_list = ["Homo sapiens","Gallus gallus"]
    result = opentree.util.get_phylo_context(species_list)
    assert result == "Tetrapods"

    species_list = ["Homo sapiens"]
    result = opentree.util.get_phylo_context(species_list)
    assert result == "Mammals"

    df = pd.DataFrame({"species":["Homo sapiens"]})
    good_species_list = [np.array(["Homo sapiens"]),
                         df.loc[:,"species"],
                         ("Homo sapiens",)]

    for g in good_species_list:
        result = opentree.util.get_phylo_context(g)
        assert result == "Mammals"

    species_list = []
    result = opentree.util.get_phylo_context(species_list)
    assert result == "All life"

    with pytest.raises(ValueError):
        species_list = ["Not really a species"]
        opentree.util.get_phylo_context(species_list)

    with pytest.raises(ValueError):
        species_list = ["Not really a species","Homo sapiens"]
        opentree.util.get_phylo_context(species_list)

    bad_species_list = [None,[1.5,2.5],1,{1:"test"},list,int,str]
    for b in bad_species_list:
        print(f"Trying bad species list {b}")
        with pytest.raises(ValueError):
            opentree.util.get_phylo_context(b)

def test_is_allowed_phylo_context():

    assert "Animals" == opentree.util.is_allowed_phylo_context("Animals")

    with pytest.raises(ValueError):
        opentree.util.is_allowed_phylo_context("not_real_context")


def test_species_to_ott():

    species_list = ["Homo sapiens",
                    "Thermus thermophilus",
                    "Methanococcus voltae",
                    "Saccharomyces cerevisiae"]

    expected_ott = [770315,276534,565131,356221]

    out = opentree.util.species_to_ott(species_list,phylo_context="All life")
    assert len(out[1]) == 0
    for i, s in enumerate(species_list):
        assert out[0][s][0] == expected_ott[i]

    bad_species_list = [None,[1.5,2.5],1,{1:"test"},list,int,str,"SPECIES"]
    for b in bad_species_list:
        print(f"Trying bad species list {b}")
        with pytest.raises(ValueError):
            opentree.util.species_to_ott(b)

    out = opentree.util.species_to_ott(["Nannospalax galili"],
                                        phylo_context="All life")
    # will throw a key error if regex not properly stripping "(in domain Bacteria)"
    out[0]["Nannospalax galili"]
    assert out[0]["Nannospalax galili"][1] == "Nannospalax galili"

    # Make sure phylo_context is really working. This species name is ambiguous;
    # this checks to see if got correct one
    out = opentree.util.species_to_ott(["Nannospalax galili"],
                                        phylo_context="Bacteria")
    assert out[0]["Nannospalax galili"][0] == 5909124
