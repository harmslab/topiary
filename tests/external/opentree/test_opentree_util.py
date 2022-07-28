import pytest
import topiary
from topiary.external.opentree.util import _validate_ott_vs_species
from topiary.external.opentree.util import ott_mrca, ott_resolvable
from topiary.external.opentree.util import species_to_ott, ott_species_tree

import ete3

import pandas as pd
import numpy as np

import re

def test__validate_ott_vs_species():

    # Nothing in
    with pytest.raises(ValueError):
        _validate_ott_vs_species(None,None)

    # Both in
    with pytest.raises(ValueError):
        _validate_ott_vs_species([],[])

    # Bad ott
    with pytest.raises(ValueError):
        _validate_ott_vs_species(ott_list=["Homo sapiens","Gallus gallus"])

    # Bad species
    with pytest.raises(ValueError):
        _validate_ott_vs_species(species_list=[111,222])

    # Good ott
    ott_list = _validate_ott_vs_species(ott_list=[111,222])
    assert len(ott_list) == 2
    assert ott_list[0] == 111
    assert ott_list[1] == 222

    # Good species
    ott_list = _validate_ott_vs_species(species_list=["Homo sapiens"])
    assert len(ott_list) == 1
    assert ott_list[0] == 770315

    # empty ott
    ott_list = _validate_ott_vs_species(ott_list=[])
    assert len(ott_list) == 0

    # empty species
    ott_list = _validate_ott_vs_species(species_list=[])
    assert len(ott_list) == 0




def test_species_to_ott():

    input_species_list = ["Homo sapiens",
                          "Thermus thermophilus",
                          "Methanococcus voltae",
                          "Saccharomyces cerevisiae"]

    expected_ott = [770315,276534,565131,356221]
    expected_taxid = [9606,274,2188,4932]

    ott_list, species_list, results = species_to_ott(input_species_list)
    assert len(ott_list) == 4
    assert len(species_list) == 4
    assert len(results) == 4
    assert np.array_equal(ott_list,expected_ott)
    assert np.array_equal(species_list,input_species_list)
    assert issubclass(type(results),dict)
    for i, s in enumerate(input_species_list):
        assert results[s]["num_matches"] == 1
        assert results[s]["taxid"] == expected_taxid[i]

    # Make sure it checks for inputs correctly
    bad_inputs = [None,[1.5,2.5],1,{1:"test"},list,int,str,"SPECIES"]
    for b in bad_inputs:
        print("Trying",b)
        with pytest.raises(ValueError):
            species_to_ott(b)

    # Make sure it handles an empty input list gracefully
    ott_list, species_list, results = species_to_ott([])
    assert len(ott_list) == 0
    assert len(species_list) == 0

    ott_list, species_list, results = species_to_ott(["Not really a species"])
    assert len(ott_list) == 1
    assert len(species_list) == 1
    assert ott_list[0] is None
    assert species_list[0] == "Not really a species"


def test_ott_species_tree():

    input_ott = [770315,276534,565131,356221]

    T, results = ott_species_tree(input_ott)
    assert issubclass(type(T),ete3.Tree)
    assert issubclass(type(results),dict)
    seen = list(results["resolved"])
    sent_in = input_ott[:]
    seen.sort()
    sent_in.sort()

    assert np.array_equal(seen,sent_in)
    assert len(results["not_resolved"]) == 0
    assert len(results["unknown_ids"]) == 0
    assert len(results["not_monophyletic"]) == 0

    off_tree = []
    for n in T.get_leaves():
        off_tree.append(int(n.name[3:]))
    off_tree.sort()
    assert np.array_equal(off_tree,sent_in)

     # Send in a bad ott
    input_ott = [770315,276534,565131,356221,9999999999999999999999]
    T, results = ott_species_tree(input_ott)

    seen = list(results["resolved"])
    sent_in = input_ott[:]
    seen.sort()
    sent_in.sort()
    assert np.array_equal(seen,sent_in[:-1])

    off_tree = []
    for n in T.get_leaves():
        off_tree.append(int(n.name[3:]))
    off_tree.sort()
    assert np.array_equal(off_tree,sent_in[:-1])

    assert len(results["not_resolved"]) == 1
    assert results["not_resolved"][0] == 9999999999999999999999
    assert len(results["unknown_ids"]) == 1
    assert results["unknown_ids"][0] == 9999999999999999999999
    assert len(results["not_monophyletic"]) == 0

    # Make sure it handles empty list gracefully
    T, results = ott_species_tree([])
    assert T is None
    assert len(results["resolved"]) == 0
    assert len(results["unknown_ids"]) == 0
    assert len(results["not_monophyletic"]) == 0
    assert len(results["not_resolved"]) == 0

     # Send in only one ott
    input_ott = [770315]
    T, results = ott_species_tree(input_ott)

    seen = list(results["resolved"])
    sent_in = input_ott[:]
    seen.sort()
    sent_in.sort()
    assert np.array_equal(seen,sent_in)

    off_tree = []
    for n in T.get_leaves():
        off_tree.append(int(n.name[3:]))
    off_tree.sort()
    assert np.array_equal(off_tree,sent_in)

    assert len(results["not_resolved"]) == 0
    assert len(results["unknown_ids"]) == 0
    assert len(results["not_monophyletic"]) == 0

     # Send in only one ott, but it's bad
    input_ott = [9999999999999999999999]
    T, results = ott_species_tree(input_ott)
    assert T is None
    assert len(results["resolved"]) == 0
    assert len(results["not_resolved"]) == 1
    assert results["not_resolved"][0] == 9999999999999999999999
    assert len(results["unknown_ids"]) == 1
    assert results["unknown_ids"][0] == 9999999999999999999999
    assert len(results["not_monophyletic"]) == 0

def test_ott_resolvable():

    some_good_ott = [770315,276534,565131,356221]

    # hybrid that is not on synthetic tree
    bad_ott = [4942641]

    resolvable = ott_resolvable(some_good_ott)
    assert np.array_equal(np.ones(len(some_good_ott),dtype=bool),resolvable)

    resolvable = ott_resolvable(bad_ott)
    assert np.array_equal(np.zeros(len(bad_ott),dtype=bool),resolvable)

    both_together = some_good_ott[:]
    both_together.extend(bad_ott)
    resolvable = ott_resolvable(both_together)
    expected = np.array([1,1,1,1,0],dtype=bool)
    assert np.array_equal(resolvable,expected)

    replicated = [770315,770315,770315,770315,770315,4942641,4942641,4942641,4942641]
    expected = np.array([1,1,1,1,1,0,0,0,0],dtype=bool)
    resolvable = ott_resolvable(replicated)
    assert np.array_equal(expected,resolvable)

    # Make sure it handes one good taxa well
    resolvable = ott_resolvable([770315])
    assert len(resolvable) == 1
    assert resolvable[0] == True

    # Make sure it handes one bad taxa well
    resolvable = ott_resolvable([4942641])
    assert len(resolvable) == 1
    assert resolvable[0] == False

    # Make sure it handles single empty list
    resolvable = ott_resolvable([])
    assert len(resolvable) == 0

    # Make sure it checks for inputs correctly
    bad_inputs = [None,1,list,int,str,"SPECIES"]
    for b in bad_inputs:
        print("Trying",b)
        with pytest.raises(ValueError):
            ott_resolvable(b)


def test_ott_mrca():

    vert = [770315, 153563, 1005914]           #["Homo sapiens","Gallus gallus","Danio rerio"]
    amniotes = [770315, 153563]                #["Homo sapiens","Gallus gallus"]
    tetrapods = [770315, 153563, 465096]       #["Homo sapiens","Gallus gallus","Xenopus laevis"]
    placental_mammals = [770315, 542509]       #["Homo sapiens","Mus musculus"]
    therian_mammals = [770315, 542509, 122362] #["Homo sapiens","Mus musculus","Monodelphis domestica"]
    plants = [309263,1075713]                  #["Arabidopsis thaliana","Pseudotsuga menziesii"]
    eukaryotes = [770315, 356221, 1075713]     #["Homo sapiens","Saccharomyces cerevisiae","Pseudotsuga menziesii"]
    human_yeast = [770315, 356221]             #["Homo sapiens","Saccharomyces cerevisiae"]
    archaea = [7000614,996421]                 #["Mancarchaeum acidiphilum","Mendosicutes"]
    not_bacteria = [770315,7000614,996421]     #["Homo sapiens","Mancarchaeum acidiphilum","Mendosicutes"]
    all_life = [770315, 474506, 7000614]       #["Homo sapiens","Escherichia coli","Mancarchaeum acidiphilum"]
    bacteria = [1090496, 474506]               #["Staphylococcus aureus","Escherichia coli"]

    # Get with default args
    out = ott_mrca(vert)
    assert out["ott_name"] == 'Euteleostomi'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 117571

    # Move up four ranks from Euteleostomi
    out = ott_mrca(vert,move_up_by=4)
    assert out["ott_name"] == 'Craniata'
    assert out["ott_rank"] == 'subphylum'
    assert out["taxid"] == 89593

    # Move waaaaay up but do not allow all life (end up on Eukaryota)
    out = ott_mrca(vert,move_up_by=10000)
    assert out["ott_name"] == 'Eukaryota'
    assert out["ott_rank"] == 'domain'
    assert out["taxid"] == 2759

    # Move waaaaay up, allowing all life
    out = ott_mrca(vert,avoid_all_life=False,move_up_by=10000)
    assert out["ott_name"] == 'life'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 1

    # Try to avoid all life, but we can't because inputs are from all three
    # domains.
    out = ott_mrca(all_life,avoid_all_life=True)
    assert out["ott_name"] == 'cellular organisms'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 131567

    # Spot check various taxonomic groups
    out = ott_mrca(amniotes)
    assert out["ott_name"] == 'Amniota'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 32524

    out = ott_mrca(tetrapods)
    assert out["ott_name"] == 'Tetrapoda'
    assert out["ott_rank"] == 'superclass'
    assert out["taxid"] == 32523

    out = ott_mrca(placental_mammals)
    assert out["ott_name"] == 'Euarchontoglires'
    assert out["ott_rank"] == 'superorder'
    assert out["taxid"] == 314146

    out = ott_mrca(therian_mammals)
    assert out["ott_name"] == 'Theria'
    assert out["ott_rank"] == 'subclass'
    assert out["taxid"] == 32525

    out = ott_mrca(human_yeast)
    assert out["ott_name"] == 'Opisthokonta'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 33154

    out = ott_mrca(plants)
    assert out["ott_name"] == 'Spermatophyta'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 58024

    out = ott_mrca(eukaryotes)
    assert out["ott_name"] == 'Eukaryota'
    assert out["ott_rank"] == 'domain'
    assert out["taxid"] == 2759

    out = ott_mrca(archaea)
    assert out["ott_name"] == 'Archaea'
    assert out["ott_rank"] == 'domain'
    assert out["taxid"] == 2157

    out = ott_mrca(bacteria)
    assert out["ott_name"] == 'Bacteria'
    assert out["ott_rank"] == 'domain'
    assert out["taxid"] == 2

    out = ott_mrca(not_bacteria)
    assert out["ott_name"] == 'cellular organisms'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 131567

    out = ott_mrca(all_life)
    assert out["ott_name"] == 'cellular organisms'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 131567

    # Make sure it checks for inputs correctly
    bad_inputs = [None,1,list,int,str,"SPECIES"]
    for b in bad_inputs:
        print("Trying",b)
        with pytest.raises(ValueError):
            ott_mrca(b)

    # Pass in bad values for move_up_by and avoid_all_life
    with pytest.raises(ValueError):
        ott_mrca(vert,move_up_by=-1)

    with pytest.raises(ValueError):
        ott_mrca(vert,move_up_by="stupid")

    with pytest.raises(ValueError):
        ott_mrca(vert,avoid_all_life="True")

    # Send in empty list
    out = ott_mrca([])
    assert out["ott_name"] == 'life'
    assert out["ott_rank"] == 'no rank'
    assert out["taxid"] == 1

    # Send in only bad ott
    with pytest.raises(ValueError):
        out = ott_mrca([99999999999999])

    # Send in placentl mammals with a bad ott
    to_test = placental_mammals[:]
    to_test.append(99999999999999)
    out = ott_mrca(to_test)
    assert out["ott_name"] == 'Euarchontoglires'
    assert out["ott_rank"] == 'superorder'
    assert out["taxid"] == 314146
