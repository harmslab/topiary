import pytest

import topiary

import numpy as np
import pandas as pd

import os, subprocess, sys

def test_main(test_dataframes,tmpdir):

    df = test_dataframes["good-df"].copy()

    # Create dataframe file
    df_file = os.path.abspath(os.path.join(tmpdir,"test-df.csv"))
    topiary.write_dataframe(df,df_file)

    # Create fasta file
    fasta_file = os.path.abspath(os.path.join(tmpdir,"some-fasta-file.fasta"))
    f = open(fasta_file,'w')
    for i in range(len(df)):
        f.write(f">{df.loc[df.index[i],'uid']}\n")
        f.write(f"{df.loc[df.index[i],'sequence']}\n")
    f.close()

    out_file = os.path.abspath(os.path.join(tmpdir,"output-df.csv"))

    # Get location of binary
    location = os.path.dirname(os.path.realpath(__file__))
    test_bin = os.path.join(location,"..","..","bin","topiary-fasta-into-dataframe")

    if os.name == "nt":
        base_cmd = [sys.executable,test_bin]
    else:
        base_cmd = [test_bin]

    # Should run but fail because no arguments
    ret = subprocess.run(base_cmd)
    assert ret.returncode != 0

    # Should run and create new csv file with a_new_column that has the
    # sequences in it
    cmd = base_cmd[:]
    cmd.extend([df_file,fasta_file,out_file,"--load_into_column","a_new_column"])
    ret = subprocess.run(cmd)
    assert ret.returncode == 0
    assert os.path.isfile(out_file)
    out_df = topiary.read_dataframe(out_file)
    assert np.array_equal(df.sequence,out_df.a_new_column)

    # Should run and create a new csv file with alignment that has the
    # sequences in it.

    # Create fasta file with identical sequence lengths
    fasta_file = os.path.abspath(os.path.join(tmpdir,"some-fasta-file.fasta"))
    f = open(fasta_file,'w')
    for i in range(len(df)):
        f.write(f">{df.loc[df.index[i],'uid']}\n")
        f.write(f"QQQQQQ\n")
    f.close()

    os.remove(out_file)

    cmd = base_cmd[:]
    cmd.extend([df_file,fasta_file,out_file])
    ret = subprocess.run(cmd)
    assert ret.returncode == 0
    assert os.path.isfile(out_file)
    out_df = topiary.read_dataframe(out_file)
    assert np.array_equal(["QQQQQQ" for _ in range(len(out_df.alignment))],
                          out_df.alignment)


    # Should run and create a new csv file with alignment that has the
    # sequences in it, with first entry keep[0] == False

    # Create fasta file with identical sequence lengths
    fasta_file = os.path.abspath(os.path.join(tmpdir,"some-fasta-file.fasta"))
    f = open(fasta_file,'w')
    for i in range(1,len(df)):
        f.write(f">{df.loc[df.index[i],'uid']}\n")
        f.write(f"QQQQQQ\n")
    f.close()

    os.remove(out_file)

    cmd = base_cmd[:]
    cmd.extend([df_file,fasta_file,out_file])
    ret = subprocess.run(cmd)
    assert ret.returncode == 0
    assert os.path.isfile(out_file)
    out_df = topiary.read_dataframe(out_file)
    assert pd.isna(out_df.loc[out_df.index[0],"alignment"])
    assert not out_df.loc[out_df.index[0],"keep"]
    for i in range(1,len(out_df)):
        assert out_df.loc[out_df.index[i],"alignment"] == "QQQQQQ"
        assert out_df.loc[out_df.index[i],"keep"]

    # Should run and create a new csv file with alignment that has the
    # sequences in it. SHould keep -- even blank ones

    os.remove(out_file)

    cmd = base_cmd[:]
    cmd.extend([df_file,fasta_file,out_file,"--load_into_column",
                "a_new_column","--unkeep_missing"])
    ret = subprocess.run(cmd)
    assert ret.returncode == 0
    assert os.path.isfile(out_file)
    out_df = topiary.read_dataframe(out_file)
    assert pd.isna(out_df.loc[out_df.index[0],"a_new_column"])
    for i in range(1,len(out_df)):
        assert out_df.loc[out_df.index[i],"a_new_column"] == "QQQQQQ"
    assert np.array_equal(out_df["keep"],np.ones(len(out_df["keep"])))
