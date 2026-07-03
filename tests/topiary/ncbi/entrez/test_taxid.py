
import pytest

import topiary
from topiary.ncbi.entrez.taxid import get_taxid

import numpy as np
import pandas as pd

@pytest.mark.run_ncbi_server
def test_get_taxid():

    expt_taxid = set(['9606', '10090'])
    obs_taxid = get_taxid(["Homo sapiens","Mus musculus"])
    assert len(set(obs_taxid).difference(expt_taxid)) == 0

    obs_taxid = get_taxid(["Homo sapiens"])
    assert obs_taxid[0] == "9606"

    obs_taxid = get_taxid("Homo sapiens")
    assert obs_taxid == "9606"

    obs_taxid = get_taxid([])
    assert len(obs_taxid) == 0

    with pytest.raises(RuntimeError):
        get_taxid(["Not a species"])

@pytest.mark.run_ncbi_server
def test_get_taxid_strains():
    # Escherichia coli (strain K12) should map to 562 (Escherichia coli) 
    # after stripping.
    obs_taxid = get_taxid("Escherichia coli (strain K12)")
    assert obs_taxid == "562"

    # Rhizobium meliloti (strain 1021) -> Rhizobium meliloti
    obs_taxid = get_taxid("Rhizobium meliloti (strain 1021)")
    assert obs_taxid == "382"


    bad_species = [None,1.5,{1:"test"},float,str]
    for b in bad_species:
        print(f"passing bad species {b}")
        with pytest.raises(ValueError):
            get_taxid(b)
