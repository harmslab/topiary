

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
