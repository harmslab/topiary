import pytest
from topiary.reports.cards.input import create_input_card
from topiary._private.supervisor import Supervisor
import os
import pandas as pd

def test_create_input_card(tiny_phylo, mocker):
    
    # Initialize supervisor from tiny_phylo dir
    calc_dir = tiny_phylo["03_ancestors"]
    sv = Supervisor(calc_dir)

    # Mock ott_to_mrca to avoid network calls
    mock_mrca = mocker.patch("topiary.opentree.ott_to_mrca", return_value={"ott_name": "TestTaxa"})

    # Test with default p_column (will use 'name' as recip_paralog is missing)
    html = create_input_card(sv)
    assert "Input" in html
    assert "Number of sequences" in html
    assert "TestTaxa" in html
    assert mock_mrca.called

    # Test with recip_paralog present
    sv_recip = Supervisor(calc_dir)
    sv_recip.df["recip_paralog"] = "paralog"
    html = create_input_card(sv_recip)
    assert "Input" in html

    # Test with explicit p_column
    html = create_input_card(sv, p_column="species")
    assert "Input" in html
    assert "Number of sequences" in html

    # Test with Null ott
    sv_no_ott = Supervisor(calc_dir)
    sv_no_ott.df["ott"] = None
    html = create_input_card(sv_no_ott)
    assert "Input" in html
    assert "Taxonomic distribution" in html

    # Test with ott column missing entirely (e.g. dataframe without species
    # information -- gene/species reconciliation was never run)
    sv_missing_ott = Supervisor(calc_dir)
    sv_missing_ott._df = sv_missing_ott.df.drop(columns=["ott"])
    html = create_input_card(sv_missing_ott)
    assert "Input" in html
    assert "Taxonomic distribution" in html