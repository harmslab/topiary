
import pytest

import topiary
from topiary.external.ncbi.entrez._get_taxid import get_taxid

import numpy as np
import pandas as pd

def test_get_taxid():

    expt_taxid = set(['10090', '9606', '7955'])
    obs_taxid = get_taxid(["Homo sapiens","Mus musculus","Danio rerio"])
    assert len(set(obs_taxid).difference(expt_taxid)) == 0

    obs_taxid = get_taxid(["Homo sapiens"])
    assert obs_taxid[0] == "9606"

    obs_taxid = get_taxid("Homo sapiens")
    assert obs_taxid == "9606"

    obs_taxid = get_taxid([])
    assert len(obs_taxid) == 0

    with pytest.raises(RuntimeError):
        get_taxid(["Not a species"])


    bad_species = [None,1.5,{1:"test"},float,str]
    for b in bad_species:
        print(f"passing bad species {b}")
        with pytest.raises(ValueError):
            get_taxid(b)
