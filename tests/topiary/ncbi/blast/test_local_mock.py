
import pytest
import os
import pandas as pd
import numpy as np
from topiary.ncbi.blast.local import local_blast, _prepare_for_blast, _construct_args, _local_blast_thread_function, _combine_hits

def test__prepare_for_blast(tmpdir, mocker):
    
    # 1. db not found
    with pytest.raises(FileNotFoundError):
        _prepare_for_blast("SEQ", "missing_db", "blastp", 100, 0.001, (11,1), 1, {})

    # Create dummy psq
    os.chdir(tmpdir)
    with open("exists.psq", "w") as f:
        f.write("data")
    
    # 2. blast_program not recognized
    with pytest.raises(ValueError):
        _prepare_for_blast("SEQ", "exists", "NOT_A_PROGRAM", 100, 0.001, (11,1), 1, {})

    # 3. blast_program not in path
    mocker.patch("subprocess.run", side_effect=FileNotFoundError)
    with pytest.raises(FileNotFoundError):
        _prepare_for_blast("SEQ", "exists", "blastp", 100, 0.001, (11,1), 1, {})
    mocker.patch("subprocess.run", return_value=None)

    # 4. Success
    seq_list, blast_args, singleton = _prepare_for_blast("SEQ", "exists", "blastp", 100, 0.001, (11,1), 1, {"task":"blastp-fast"})
    assert seq_list == ["SEQ"]
    assert "blastp" in blast_args
    assert "-db" in blast_args
    assert "exists" in blast_args
    assert singleton is True

def test__construct_args():
    
    # manual_num_cores is fixed below so this test does not depend on the
    # number of cores actually available on the machine running it (this
    # varies across CI runners and caused flaky failures on macOS, which
    # historically has fewer cores on GitHub-hosted runners than Linux).

    # Hit line 166: max_useful_threads > num_threads
    seq_list_10 = ["S"] * 10
    blast_args = ["blastp", "-db", "db"]
    kwargs_list, num_threads = _construct_args(seq_list_10, blast_args, block_size=2, num_threads=2, manual_num_cores=8)
    assert num_threads == 2

    # 3 sequences, block_size 2 -> 1 window, remainder 1. 1/2 = 0.5 >= 0.5 -> TRUE.
    # HOWEVER: since num_threads=4 > 1 (max_useful_threads), it will re-adjust
    # block_size to 1 to try to use more threads.
    seq_list_3 = ["S"] * 3
    kwargs_list, num_threads = _construct_args(seq_list_3, blast_args, block_size=2, num_threads=4, manual_num_cores=8)
    assert len(kwargs_list) == 3

    # Hit line 192: counter >= len(windows)
    # 12 sequences, block_size 5 -> num_windows 2, remainder 2.
    # 2/5 = 0.4 < 0.5 -> else block.
    # windows = [5, 5].
    # Loop 1: windows[0]=6, rem=1, counter=1.
    # Loop 2: windows[1]=6, rem=0, counter=2. counter >= 2 -> counter=0.
    seq_list_12 = ["S"] * 12
    kwargs_list, num_threads = _construct_args(seq_list_12, blast_args, block_size=5, num_threads=4, manual_num_cores=8)
    assert len(kwargs_list) == 4
    assert kwargs_list[0]["index"] == (0, 3)
    assert kwargs_list[1]["index"] == (3, 6)
    assert kwargs_list[2]["index"] == (6, 9)
    assert kwargs_list[3]["index"] == (9, 12)

def test__local_blast_thread_function(tmpdir, mocker):
    
    os.chdir(tmpdir)
    mocker.patch("topiary._private.interface.launch")
    mocker.patch("os.remove")
    
    # Mock read_blast_xml
    df = pd.DataFrame({"query":["count0"], "accession":["ACC1"]})
    mocker.patch("topiary.ncbi.blast.local.read_blast_xml", return_value=([df], ["fake.xml"]))
    
    seq_list = ["ACGT"]
    index = (0, 1)
    blast_args = ["blastp", "-db", "db"]
    
    # 1. Success
    out_df = _local_blast_thread_function(seq_list, index, blast_args, keep_blast_xml=False)
    assert out_df.equals(df)
    
    # 2. FileNotFoundError from internal blast
    mocker.patch("topiary.ncbi.blast.local.read_blast_xml", side_effect=FileNotFoundError)
    with pytest.raises(RuntimeError):
        _local_blast_thread_function(seq_list, index, blast_args, keep_blast_xml=False)

def test__combine_hits():
    
    df1 = pd.DataFrame({"query":["count0", "count1"], "accession":["ACC1", "ACC2"]})
    df2 = pd.DataFrame({"query":["count2"], "accession":[np.nan]}) # No hits
    
    # 1. Multiple hits
    combined = _combine_hits([df1, df2], return_singleton=False)
    assert len(combined) == 3
    
    # 2. Singleton
    combined = _combine_hits([df1.iloc[[0]]], return_singleton=True)
    assert isinstance(combined, pd.DataFrame)

def test_local_blast(tmpdir, mocker):
    
    os.chdir(tmpdir)
    with open("exists.psq", "w") as f:
        f.write("data")
        
    mocker.patch("subprocess.run", return_value=None)
    
    # Mock threading.thread_manager to return our dummy hits
    df = pd.DataFrame({"query":["count0"], "accession":["ACC1"]})
    mocker.patch("topiary._private.threads.thread_manager", return_value=[df])
    
    # 1. Successful run
    result = local_blast("ACGT", "exists")
    assert isinstance(result, pd.DataFrame)
    assert result["accession"].iloc[0] == "ACC1"
