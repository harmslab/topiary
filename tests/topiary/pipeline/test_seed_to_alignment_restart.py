import pytest
import topiary
import os
import shutil
import pandas as pd
import glob
from topiary.pipeline.seed_to_alignment import seed_to_alignment

def test_seed_to_alignment_restart_with_xml(tmpdir, mocker):
    """
    Test that seed_to_alignment correctly picks up existing XML files on restart.
    """
    
    # Mock necessary functions to avoid actual BLAST/network/filesystem calls
    mocker.patch("topiary.pipeline.seed_to_alignment.installed.validate_stack")
    mocker.patch("topiary.pipeline.seed_to_alignment.check.check_bool", side_effect=lambda x, y: x)
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.io.read_seed", 
                 return_value=(pd.DataFrame({"species":["spec1"], "name":["name1"], "aliases":["alias1"], "sequence":["AAAA"]}), 
                               ["spec1"], ["para*"], True))
    
    mock_df_from_seed = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.df_from_seed", 
                                     return_value=(pd.DataFrame(), ["spec1"], ["para*"], True))
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.write_dataframe")
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.read_dataframe", return_value=pd.DataFrame())
    
    # Recip blast mock
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.recip_blast", return_value=pd.DataFrame())
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.ncbi.get_taxid", return_value="12345")
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.ncbi.get_proteome", return_value="proteome.fasta")
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.ncbi.make_blast_db")
    
    # Quality/Align mocks
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.quality.shrink_dataset", return_value=pd.DataFrame())
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.muscle.align", return_value=pd.DataFrame())
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.quality.polish_alignment", return_value=pd.DataFrame())
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.write_fasta")

    # Set up the test directory
    out_dir = str(tmpdir.mkdir("test_out"))
    seed_file = str(tmpdir.join("seed.csv"))
    with open(seed_file, "w") as f:
        f.write("test")

    # Create a dummy XML file in the output directory
    xml_prefix = "01_topiary-initial-blast"
    xml_file = os.path.join(out_dir, f"{xml_prefix}_RANDOM_ncbi-blast-result.xml")
    with open(xml_file, "w") as f:
        f.write("<blast>dummy</blast>")

    # Run seed_to_alignment with restart
    # We need to mock os.path.exists to return True for the directory and the xml, but False for the CSV
    original_exists = os.path.exists
    def side_effect_exists(path):
        path = str(path)
        if "03_initial-dataframe.csv" in path:
            return False
        if "02_downloaded-sequences.csv" in path:
            return False
        if path == out_dir:
            return True
        if "seed.csv" in path:
            return True
        if xml_prefix in path:
            return True
        return original_exists(path)
    
    mocker.patch("os.path.exists", side_effect=side_effect_exists)
    mocker.patch("os.path.isdir", return_value=True)
    mocker.patch("os.chdir")
    mocker.patch("os.mkdir")
    mocker.patch("shutil.copy")

    def side_effect_glob(pattern):
        if xml_prefix in pattern:
            return [os.path.basename(xml_file)]
        return []
    mocker.patch("glob.glob", side_effect=side_effect_glob)

    # Run the function
    seed_to_alignment(seed_file, out_dir=out_dir, restart=True, local_blast_db="some_db")

    # Check that df_from_seed was called with blast_xml containing our dummy file
    # and local_blast_db being None
    args, kwargs = mock_df_from_seed.call_args
    assert kwargs["blast_xml"] == [os.path.basename(xml_file)]
    assert kwargs["local_blast_db"] is None
    assert kwargs["ncbi_blast_db"] is None
    assert kwargs["intermediate_file"] == "02_downloaded-sequences.csv"
    assert kwargs["name_prefix"] == xml_prefix

def test_seed_to_alignment_restart_with_intermediate_df(tmpdir, mocker):
    """
    Test that seed_to_alignment correctly picks up intermediate downloaded sequences on restart.
    """
    
    mocker.patch("topiary.pipeline.seed_to_alignment.installed.validate_stack")
    mocker.patch("topiary.pipeline.seed_to_alignment.check.check_bool", side_effect=lambda x, y: x)
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.io.read_seed", 
                 return_value=(pd.DataFrame({"species":["spec1"], "name":["name1"], "aliases":["alias1"], "sequence":["AAAA"]}), 
                               ["spec1"], ["para*"], True))
    
    mock_df_from_seed = mocker.patch("topiary.pipeline.seed_to_alignment.topiary.df_from_seed", 
                                     return_value=(pd.DataFrame(), ["spec1"], ["para*"], True))
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.write_dataframe")
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.read_dataframe", return_value=pd.DataFrame())
    
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.recip_blast", return_value=pd.DataFrame())
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.ncbi.get_taxid", return_value="12345")
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.ncbi.get_proteome", return_value="proteome.fasta")
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.ncbi.make_blast_db")
    
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.quality.shrink_dataset", return_value=pd.DataFrame())
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.muscle.align", return_value=pd.DataFrame())
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.quality.polish_alignment", return_value=pd.DataFrame())
    mocker.patch("topiary.pipeline.seed_to_alignment.topiary.write_fasta")

    out_dir = str(tmpdir.mkdir("test_out"))
    seed_file = str(tmpdir.join("seed.csv"))
    with open(seed_file, "w") as f:
        f.write("test")

    # Create an intermediate CSV file
    intermediate_file = os.path.join(out_dir, "02_downloaded-sequences.csv")
    with open(intermediate_file, "w") as f:
        f.write("accession,sequence\nACC1,AAAA")

    def side_effect_exists(x):
        x = str(x)
        if "03_initial-dataframe.csv" in x:
            return False
        if "02_downloaded-sequences.csv" in x:
            return True
        return True

    mocker.patch("os.path.exists", side_effect=side_effect_exists)
    mocker.patch("os.path.isdir", return_value=True)
    mocker.patch("os.chdir")
    mocker.patch("os.mkdir")
    mocker.patch("shutil.copy")

    seed_to_alignment(seed_file, out_dir=out_dir, restart=True, local_blast_db="some_db")

    args, kwargs = mock_df_from_seed.call_args
    # On this restart, blast_xml should be None because we found the intermediate CSV
    assert kwargs["blast_xml"] is None
    assert kwargs["local_blast_db"] is None
    assert kwargs["ncbi_blast_db"] is None
    assert kwargs["intermediate_file"] == "02_downloaded-sequences.csv"
