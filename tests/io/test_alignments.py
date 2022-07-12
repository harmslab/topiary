import pytest

from topiary import io
import numpy as np
import pandas as pd

import os, shutil, re


def _validate_seq_writer():
    """
    This function is tested implicitly by test_write_fasta and test_write_phy.
    I'm using those unit tests because I can validate written output.
    """

    pass

def test_write_fasta(test_dataframes,tmpdir):

    def _check_output_file(out_file,num_columns,check_length=None):

        f = open(out_file)
        lines = f.readlines()
        f.close()

        has_gap = False
        for l in lines:
            if l.startswith(">"):
                col = l[1:].split("|")
                assert len(col) == num_columns
                assert len(col[0]) == 10
            else:
                if re.search("-",l):
                    has_gap = True

        if check_length is not None:
            assert len(lines) == check_length*2

        return has_gap

    df = test_dataframes["good-df_with-good-alignment"]

    # Basic read/write
    out_file = os.path.join(tmpdir,"output.fasta")
    io.write_fasta(df,out_file)
    _check_output_file(out_file,num_columns=3)

    # -------------------------------------------------------------------------
    # output files

    bad_output_files = [None,1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                        str,(1,2)]
    for b in bad_output_files:
        with pytest.raises(ValueError):
            io.write_fasta(df,b)

    # -------------------------------------------------------------------------
    # seq_column

    os.remove(out_file)
    io.write_fasta(df,out_file=out_file)
    has_gaps = _check_output_file(out_file,num_columns=3)
    assert not has_gaps

    os.remove(out_file)
    io.write_fasta(df,out_file=out_file,seq_column="alignment")
    has_gaps = _check_output_file(out_file,num_columns=3)
    assert has_gaps

    bad_seq_columns = [None,1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                        str,(1,2),"not_in_df"]
    for b in bad_seq_columns:
        with pytest.raises(ValueError):
            io.write_fasta(df,out_file=out_file,seq_column=b)

    # -------------------------------------------------------------------------
    # columns

    bad_columns = [None,1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2)]
    for b in bad_columns:
        with pytest.raises(ValueError):
            io.write_fasta(df,out_file=out_file,label_columns=b)

    os.remove(out_file)
    io.write_fasta(df,out_file,label_columns=["start"])
    has_gaps = _check_output_file(out_file,num_columns=2)

    # -------------------------------------------------------------------------
    # write_only_keepers

    bad_write_only = [None,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2),"something"]
    for b in bad_write_only:
        with pytest.raises(ValueError):
            io.write_fasta(df,out_file=out_file,write_only_keepers=b)

    no_keep_df = df.copy()
    no_keep_df.loc[1:,"keep"] = False
    os.remove(out_file)
    io.write_fasta(no_keep_df,out_file,write_only_keepers=True)
    has_gaps = _check_output_file(out_file,num_columns=3,check_length=1)

    no_keep_df.keep = False
    os.remove(out_file)
    with pytest.raises(RuntimeError):
        io.write_fasta(no_keep_df,out_file,write_only_keepers=True)

    # -------------------------------------------------------------------------
    # empty_char

    bad_empty = [1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                 str,(1,2),True]
    for b in bad_empty:
        with pytest.raises(ValueError):
            io.write_fasta(df,out_file=out_file,empty_char=b)

    # Should throw runtime error because sequences all look empty!
    with pytest.raises(RuntimeError):
        io.write_fasta(df,out_file=out_file,empty_char="ACDEFGHIKLMNPQRSTVWYZ-")

    # -------------------------------------------------------------------------
    # clean_sequence

    bad_clean = [None,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2),"something"]
    for b in bad_clean:
        with pytest.raises(ValueError):
            io.write_fasta(df,out_file=out_file,clean_sequence=b)

    to_clean_df = df.copy()
    to_clean_df.loc[:,"sequence"] = "STUPIDX?STUPID"
    to_clean_df.loc[:,"length"] = len("STUPIDX?STUPID")
    io.write_fasta(to_clean_df,out_file=out_file,clean_sequence=True)
    f = open(out_file)
    lines = f.readlines()
    f.close()
    assert lines[1].strip() == "ST-PID--ST-PID"

    # -------------------------------------------------------------------------
    # overwrite

    bad_overwrite = [None,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2),"something"]
    for b in bad_overwrite:
        with pytest.raises(ValueError):
            io.write_fasta(df,out_file=out_file,overwrite=b)

    os.remove(out_file)
    io.write_fasta(df,out_file)
    # Should throw error because overwrite is False by default
    with pytest.raises(FileExistsError):
        io.write_fasta(df,out_file)

    # Should throw error
    with pytest.raises(FileExistsError):
        io.write_fasta(df,out_file,overwrite=False)

    # Should work because we have the file in place
    io.write_fasta(df,out_file,overwrite=True)

def test_write_phy(test_dataframes,tmpdir):

    def _check_output_file(out_file,check_length=None):

        f = open(out_file)
        lines = f.readlines()
        f.close()

        header = lines[0].split()
        num_seqs = int(header[0])
        seq_length = int(header[1])

        if check_length is not None:
            assert num_seqs == check_length

        counter = 0
        has_gap = False
        for l in lines[2:]:

            if counter == 0:
                assert len(l.strip()) == 10
                counter += 1
            else:
                assert len(l.strip()) == seq_length
                if re.search("-",l):
                    has_gap = True

                counter = 0

        return has_gap

    df = test_dataframes["good-df_with-good-alignment"]

    # Basic read/write
    out_file = os.path.join(tmpdir,"output.fasta")
    io.write_phy(df,out_file)
    _check_output_file(out_file)

    # -------------------------------------------------------------------------
    # output files

    bad_output_files = [None,1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                        str,(1,2)]
    for b in bad_output_files:
        with pytest.raises(ValueError):
            io.write_phy(df,b)

    # -------------------------------------------------------------------------
    # seq_column

    os.remove(out_file)
    io.write_phy(df,out_file=out_file,seq_column="sequence")
    has_gaps = _check_output_file(out_file)
    assert not has_gaps

    os.remove(out_file)
    io.write_phy(df,out_file=out_file,seq_column="alignment")
    has_gaps = _check_output_file(out_file)
    assert has_gaps

    bad_seq_columns = [None,1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                        str,(1,2),"not_in_df"]
    for b in bad_seq_columns:
        with pytest.raises(ValueError):
            io.write_phy(df,out_file=out_file,seq_column=b)

    # Send in alignment column with a short sequence (not all same length)
    bad_align_df = df.copy()
    bad_align_df.loc[0,"alignment"] = "MTG"
    os.remove(out_file)
    with pytest.raises(ValueError):
        io.write_phy(bad_align_df,out_file=out_file,seq_column="alignment")

    # -------------------------------------------------------------------------
    # write_only_keepers

    bad_write_only = [None,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2),"something"]
    for b in bad_write_only:
        with pytest.raises(ValueError):
            io.write_phy(df,out_file=out_file,write_only_keepers=b)

    no_keep_df = df.copy()
    no_keep_df.loc[1:,"keep"] = False
    io.write_phy(no_keep_df,out_file,write_only_keepers=True)
    has_gaps = _check_output_file(out_file,check_length=1)

    no_keep_df.keep = False
    os.remove(out_file)
    with pytest.raises(ValueError):
        io.write_phy(no_keep_df,out_file,write_only_keepers=True)

    # -------------------------------------------------------------------------
    # empty_char

    bad_empty = [1.0,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                 str,(1,2),True]
    for b in bad_empty:
        with pytest.raises(ValueError):
            io.write_phy(df,out_file=out_file,empty_char=b)

    # Should throw runtime error because sequences all look empty!
    with pytest.raises(RuntimeError):
        io.write_phy(df,out_file=out_file,empty_char="ACDEFGHIKLMNPQRSTVWYZ-")

    # -------------------------------------------------------------------------
    # clean_sequence

    bad_clean = [None,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2),"something"]
    for b in bad_clean:
        with pytest.raises(ValueError):
            io.write_phy(df,out_file=out_file,clean_sequence=b)

    to_clean_df = df.copy()
    to_clean_df.loc[0,"alignment"] = "?LUFLFF---"
    io.write_phy(to_clean_df,out_file=out_file,clean_sequence=True)
    f = open(out_file)
    lines = f.readlines()
    f.close()
    assert lines[3].strip() == "-L-FLFF---"

    # -------------------------------------------------------------------------
    # overwrite

    bad_overwrite = [None,[1,3,4],{"test":1},pd.DataFrame({"test":[1]}),
                   str,(1,2),"something"]
    for b in bad_overwrite:
        with pytest.raises(ValueError):
            io.write_phy(df,out_file=out_file,overwrite=b)

    os.remove(out_file)
    io.write_phy(df,out_file)
    # Should throw error because overwrite is False by default
    with pytest.raises(FileExistsError):
        io.write_phy(df,out_file)

    # Should throw error
    with pytest.raises(FileExistsError):
        io.write_phy(df,out_file,overwrite=False)

    # Should work because we have the file in place
    io.write_phy(df,out_file,overwrite=True)


def test_read_fasta_into(test_dataframes,tmpdir):

    df = test_dataframes["good-df"].copy()

    # Make sure alignment gets read into test_column (tests load_into_column)
    fasta = os.path.join(tmpdir,"tmp.fasta")
    io.write_fasta(df,fasta)
    new_df = io.read_fasta_into(df,fasta,load_into_column="test_column")
    assert np.array_equal(new_df.loc[:,"test_column"],df.loc[:,"sequence"])

    # Load in lines of fasta file
    fasta_lines = []
    for i in range(len(df)):
        fasta_lines.append(f">{df.uid.iloc[i]}\n")
        fasta_lines.append(f"{df.sequence.iloc[i]}\n")
    new_df = io.read_fasta_into(df,fasta_lines,load_into_column="test_column")
    assert np.array_equal(new_df.loc[:,"test_column"],df.loc[:,"sequence"])

    # Load in lines of fasta file, no \n
    fasta_lines = []
    for i in range(len(df)):
        fasta_lines.append(f">{df.uid.iloc[i]}")
        fasta_lines.append(f"{df.sequence.iloc[i]}")
    new_df = io.read_fasta_into(df,fasta_lines,load_into_column="test_column")
    assert np.array_equal(new_df.loc[:,"test_column"],df.loc[:,"sequence"])

    # Pass in file without uid
    fasta = os.path.join(tmpdir,"tmp.fasta")
    f = open(fasta,'w')
    f.write(">bad\nQQQQQQQQQQ\n")
    f.close()
    with pytest.raises(ValueError):
        new_df = io.read_fasta_into(df,fasta,load_into_column="test_column")

    # Pass in file without uid, another format
    fasta = os.path.join(tmpdir,"tmp.fasta")
    f = open(fasta,'w')
    f.write(">bad|bad\nQQQQQQQQQQ\n")
    f.close()
    with pytest.raises(ValueError):
        new_df = io.read_fasta_into(df,fasta,load_into_column="test_column")

    # testing unkeep missing

    # Make sure if we don't write out a sequence, when we load the alignment
    # back in the dataframe drops the unwritten sequence
    df = test_dataframes["good-df"].copy()
    test_df = df.copy()
    test_df.loc[df.index[0],"keep"] = False
    io.write_fasta(test_df,fasta,overwrite=True)
    new_df = io.read_fasta_into(df,fasta,load_into_column="test_column")
    assert new_df.loc[new_df.index[0],"keep"] == False

    # But still keep if unkeep_missing == False
    df = test_dataframes["good-df"].copy()
    test_df = df.copy()
    test_df.loc[df.index[0],"keep"] = False
    io.write_fasta(test_df,fasta,overwrite=True)
    new_df = io.read_fasta_into(df,fasta,load_into_column="test_column",unkeep_missing=False)
    assert new_df.loc[new_df.index[0],"keep"] == True
