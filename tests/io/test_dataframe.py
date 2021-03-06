import pytest

from topiary import io
import numpy as np
import pandas as pd

import warnings, os, shutil, re

def test_read_dataframe(dataframe_good_files,test_dataframes):
    """
    Test read dataframe function.
    """

    ref_df = test_dataframes["good-df"]

    for f in dataframe_good_files:

        # Print f in case parsing dies... so we know which file causes failure.
        print(f)

        # Read file and make sure it does not throw warning.
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

    # Make sure dies with useful error
    bad_inputs = [1,-1,1.5,None,False,pd.DataFrame]
    for b in bad_inputs:
        with pytest.raises(ValueError):
            io.read_dataframe(b)

    # Make sure raises file not found if a file is not passed
    with pytest.raises(FileNotFoundError):
        io.read_dataframe("not_really_a_file.txt")

def test_write_dataframe(test_dataframes,tmpdir):

    df = test_dataframes["good-df"]

    bad_df = [pd.DataFrame,pd.DataFrame({"test":[1]}),None,1,"string",str]
    for b in bad_df:
        with pytest.raises(ValueError):
            io.write_dataframe(b,"output_file.csv")

    bad_out_file = [pd.DataFrame,pd.DataFrame({"test":[1]}),None,1,str]
    for b in bad_out_file:
        with pytest.raises(ValueError):
            io.write_dataframe(df,b)

    def _check_written_out(df,out,sep):
        assert os.path.isfile(out)
        with open(out) as f:
            for line in f:
                assert len(line.split(sep)) == len(df.columns)

    out = os.path.join(tmpdir,"stupid.csv")
    f = open(out,"w")
    f.write("stuff")
    f.close()
    with pytest.raises(FileExistsError):
        io.write_dataframe(df,out_file=out)

    io.write_dataframe(df,out_file=out,overwrite=True)
    _check_written_out(df,out,",")

    # Write out as csv with non-standard extension
    out = os.path.join(tmpdir,"some_file.txt")
    io.write_dataframe(df,out_file=out)
    _check_written_out(df,out,",")

    out = os.path.join(tmpdir,"some_file.tsv")
    io.write_dataframe(df,out_file=out)
    _check_written_out(df,out,"\t")

    out = os.path.join(tmpdir,"some_file.xlsx")
    io.write_dataframe(df,out_file=out)
    assert os.path.exists(out)
