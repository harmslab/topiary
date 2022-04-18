import pytest

from topiary import io
import numpy as np
import pandas as pd

import warnings

def test_check_topiary_dataframe(test_dataframes):
    """
    Test check for topiary dataframe.
    """

    # Make sure reads good dataframe without mangling
    good_df = test_dataframes["good-df"]
    df = io.check_topiary_dataframe(test_dataframes["good-df"])
    assert np.sum(np.asarray(good_df == df) == False) == 0

    # Check add uid column. Should warn that it's editing uid
    with pytest.warns():
        df = io.check_topiary_dataframe(test_dataframes["no-uid"])
    assert len(np.unique(df.uid)) == len(df)

    # Check replace bad uid. Should warn that it's editing uid
    with pytest.warns():
        df = io.check_topiary_dataframe(test_dataframes["bad-uid"])
    for d in df.uid:
        assert type(d) is str
        assert len(d) == 10

    # Check make uid unique. Should warn that it's editing uid
    with pytest.warns():
        df = io.check_topiary_dataframe(test_dataframes["duplicate-uid"])
    assert len(np.unique(df.uid)) == len(df)

    # Make sure it properly looks for required columns
    required_col = ["species","name","sequence"]
    for c in required_col:
        bad_df = good_df.copy()
        bad_df = bad_df.drop(columns=[c])
        with pytest.raises(ValueError):
            io.check_topiary_dataframe(bad_df)


def test_read_dataframe(dataframe_good_files,test_dataframes):
    """
    Test read dataframe function.
    """

    ref_df = test_dataframes["good-df"]

    for f in dataframe_good_files:

        # Print f in case parsing dies... so we know which file causes failure.
        print(f)

        # Read file and make sure it does not throuw warning.
        with warnings.catch_warnings():
            warnings.simplefilter("error")

            df = io.read_dataframe(f)
            assert len(df) == len(ref_df)

            # Make sure expected columns are present
            df.uid
            df.keep

    # Check reading a dataframe
    df = io.read_dataframe(ref_df)
    assert len(df) == len(ref_df)
    assert df is not ref_df # make sure it's a copy
    df.uid
    df.keep

    # Check "remove_extra_index" flag
    ref_df.to_csv("junk.csv")
    df = io.read_dataframe("junk.csv",remove_extra_index=False)
    assert df.columns[0].startswith("Unnamed:")
    assert len(df) == len(ref_df)

    df = io.read_dataframe("junk.csv",remove_extra_index=True)
    assert not df.columns[0].startswith("Unnamed:")
    assert len(df) == len(ref_df)

    # Make sure dies with useful error
    bad_inputs = [1,-1,1.5,None,False,pd.DataFrame]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            io.read_dataframe(b)

    # Make sure raises file not found if a file is not passed
    with pytest.raises(FileNotFoundError):
        io.read_dataframe("not_really_a_file.txt")



# def test_ncbi_blast_xml_to_df():
#
#     pass
#     #xml_files,
#     # aliases=None,
#     #phylo_context="All life"
#
# def test_write_fasta():
#     pass
#     #df,out_file,seq_column="sequence",seq_name="pretty",
#     #write_only_keepers=True,empty_char="X-?",clean_sequence=False)
#
# def test_write_phy():
#     pass
#     #df,out_file,seq_column="sequence",
#     #          write_only_keepers=True,
#     #          empty_char="X-?",
#     #          clean_sequence=False):
#
# def test_load_fasta():
#     pass
#     #df,fasta_file,load_into_column="alignment",empty_char="X-?",unkeep_missing=True):
