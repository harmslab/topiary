import pytest
import pandas as pd
import os
import io
from topiary.io.uniprot_to_seed import uniprot_to_seed
from topiary.cli_scripts.uniprot_to_seed import main as cli_main

def test_uniprot_to_seed(tmpdir):
    # Create a sample fasta file
    fasta_content = (
        ">sp|P03023|LACI_ECOLI Lactose operon repressor OS=Escherichia coli (strain K12) OX=83333 GN=lacI PE=1 SV=3\n"
        "MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQ\n"
        "SLLIGVATSSLALHAPSQIVAAIKSRADQLGASVVVSMVERSGVEACKAAVHNLLAQRVS\n"
        ">tr|A0A0B4J1G0|A0A0B4J1G0_9GAMM Uncharacterized protein OS=Vibrio sp. OX=671 GN=vibr PE=4 SV=1\n"
        "MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQ\n"
        ">sp|P12345|FAKE_PROT Fake protein with no GN OS=Fake species OX=123 PE=1 SV=1\n"
        "AAAAA\n"
        # Shared protein name between these two -> hits deduplication logic
        ">sp|P00001|PROT_A Common Protein Name OS=Species A OX=1 GN=geneA PE=1 SV=1\n"
        "CCCCC\n"
        ">sp|P00002|PROT_B Common Protein Name OS=Species B OX=2 GN=geneB PE=1 SV=1\n"
        "DDDDD\n"
        # Entry with no protein name -> covers line 103 (continue in alias loop)
        ">sp|P67890|NO_NAME_PROT OS=Null Species OX=0 PE=1 SV=1\n"
        "EEEEE\n"
        # Entry with no KV pairs -> covers lines 72-73
        ">sp|P55555|NO_KV_PROT Some protein name with no KV pairs\n"
        "FFFFF\n"
        # This one has no pipes -> should hit line 52
        ">malformed_header\n"
        "BBBBB\n"
    )
    
    fasta_path = os.path.join(tmpdir, "test.fasta")
    with open(fasta_path, "w") as f:
        f.write(fasta_content)
        
    # Run uniprot_to_seed
    # Note: read_with_read_seed=False because we don't want to hit OpenTree API in tests
    # if possible, or we should mock it. topiary.io.read_seed calls species_to_ott.
    
    df = uniprot_to_seed(fasta_path, read_with_read_seed=False)
    
    assert len(df) == 7
    
    # Check first entry
    assert df.loc[0, "species"] == "Escherichia coli (strain K12)"
    assert df.loc[0, "name"] == "lacI"
    assert "P03023" in df.loc[0, "aliases"]
    assert "LACI_ECOLI" in df.loc[0, "aliases"]
    assert "Lactose operon repressor" in df.loc[0, "aliases"]
    
    # Check second entry
    assert df.loc[1, "species"] == "Vibrio sp."
    assert df.loc[1, "name"] == "vibr"
    
    # Check third entry (no GN)
    assert df.loc[2, "species"] == "Fake species"
    assert df.loc[2, "name"] == "FAKE_PROT"

    # Check entries with shared protein name
    assert "Common Protein Name" in df.loc[3, "aliases"]
    assert "Common Protein Name" not in df.loc[4, "aliases"]
    
def test_cli(tmpdir):
    fasta_content = (
        ">sp|P03023|LACI_ECOLI Lactose operon repressor OS=Escherichia coli (strain K12) OX=83333 GN=lacI PE=1 SV=3\n"
        "MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQ\n"
    )
    fasta_path = os.path.join(tmpdir, "test.fasta")
    with open(fasta_path, "w") as f:
        f.write(fasta_content)
        
    out_path = os.path.join(tmpdir, "out.csv")
    
    # Simple mock for read_seed to avoid network calls
    def mock_read_seed(df, **kwargs):
        # Add basic columns if missing (topiary-specific)
        if "uid" not in df.columns:
            df["uid"] = ["abc"] * len(df)
        if "keep" not in df.columns:
            df["keep"] = [True] * len(df)
        return df, None, None, None

    from unittest.mock import patch
    
    with patch("topiary.io.uniprot_to_seed.read_seed", side_effect=mock_read_seed):
        with patch("topiary.write_dataframe") as mock_write:
            # Test with explicit argv
            cli_main(["--out", out_path, fasta_path])
            assert mock_write.called
            
            # Test with argv=None (covers lines 19-20 in cli script)
            with patch("sys.argv", ["topiary-uniprot-to-seed", "--out", out_path, fasta_path]):
                cli_main()
                assert mock_write.called

def test_uniprot_to_seed_with_read_seed(tmpdir):
    # Test that it calls read_seed correctly
    fasta_content = (
        ">sp|P03023|LACI_ECOLI Lactose operon repressor OS=Escherichia coli OX=83333 GN=lacI PE=1 SV=3\n"
        "MKPVTLYDVAEYAGVSYQTVSRVVNQASHVSAKTREKVEAAMAELNYIPNRVAQQLAGKQ\n"
    )
    fasta_path = os.path.join(tmpdir, "test.fasta")
    with open(fasta_path, "w") as f:
        f.write(fasta_content)
        
    from unittest.mock import patch
    
    # Mock read_seed to return a tuple
    mock_ret = (pd.DataFrame({"species":["Escherichia coli"]}), None, None, None)
    
    with patch("topiary.io.uniprot_to_seed.read_seed", return_value=mock_ret):
        df = uniprot_to_seed(fasta_path, read_with_read_seed=True)
        assert isinstance(df, pd.DataFrame)
