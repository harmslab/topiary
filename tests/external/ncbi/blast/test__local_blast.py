
import pytest
from conftest import get_public_param_defaults

import topiary
from topiary.external.ncbi.blast._local_blast import local_blast as local_blast
from topiary.external.ncbi.blast._local_blast import _prepare_for_blast as _pfb
from topiary.external.ncbi.blast._local_blast import _construct_args as _ca
from topiary.external.ncbi.blast._local_blast import _combine_hits

import Bio.Blast.Applications as apps

import numpy as np
import pandas as pd

import copy, os
import multiprocessing as mp


def test__prepare_blast(test_dataframes,tmpdir):

    # Make a fake blast db so code passes "file exists" check
    f = open(os.path.join(tmpdir,"GRCh38.psq"),"w")
    f.write("\n")
    f.close()
    fake_blast_db = os.path.join(tmpdir,"GRCh38")

    default_kwargs = get_public_param_defaults(local_blast,_pfb)
    default_kwargs["db"] = fake_blast_db
    default_kwargs["test_skip_blast_program_check"] = True
    default_kwargs["kwargs"] = {}
    df = test_dataframes["good-df"].copy()
    default_kwargs["sequence"] = df.sequence

    # Skip the check for a working blast program.
    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**default_kwargs)

    # Should be a list of length 1 holding a single dictionary
    assert type(sequence_list) is list
    assert len(sequence_list) == 5
    assert blast_function is apps.NcbiblastpCommandline
    assert type(blast_kwargs) is dict
    assert return_singleton == False

    # Make sure it got the sequences right
    assert np.array_equal(sequence_list,df.sequence)

    # Make sure the rest of the arguments were done properly
    assert blast_kwargs["max_target_seqs"] == 100
    assert blast_kwargs["threshold"] == 0.01
    assert blast_kwargs["gapopen"] == 11
    assert blast_kwargs["gapextend"] == 1

    # -------------------------------------------------------------------------
    # test sequence bits

    kwargs = copy.deepcopy(default_kwargs)

    kwargs["sequence"] = "test"
    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert len(sequence_list) == 1
    assert sequence_list[0] == "test"
    assert return_singleton == True

    kwargs["sequence"] = ["test"]
    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert len(sequence_list) == 1
    assert sequence_list[0] == "test"
    assert return_singleton == False

    kwargs["sequence"] = ["test","this"]
    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert len(sequence_list) == 2
    assert sequence_list[0] == "test"
    assert sequence_list[1] == "this"
    assert return_singleton == False

    # Send in bad sequences
    bad_sequences = [1,False,1.0,-1,None,str,""]
    for b in bad_sequences:
        print("passing bad sequence:",b)

        kwargs["sequence"] = b
        with pytest.raises(ValueError):
            sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)

        kwargs["sequence"] = [b]
        with pytest.raises(ValueError):
            sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)

    # -------------------------------------------------------------------------
    # test db

    kwargs = copy.deepcopy(default_kwargs)

    kwargs["db"] = "not_really_db"
    with pytest.raises(FileNotFoundError):
        sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)

    bad_db = [1,False,1.0,-1,None,str,"",{}]
    for b in bad_db:
        print("passing bad db:",b)
        kwargs["db"] = b
        with pytest.raises(FileNotFoundError):
            sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)

    # -------------------------------------------------------------------------
    # blast_program

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_function == apps.NcbiblastpCommandline

    kwargs["blast_program"] = "tblastx"
    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_function == apps.NcbitblastxCommandline

    bad_pg = ["not_really_program",1,False,1.0,-1,None,str,"",{}]
    for b in bad_pg:
        print("passing bad blast_program:",b)
        kwargs["blast_program"] = b
        with pytest.raises(ValueError):
            sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)

    # -------------------------------------------------------------------------
    # histlist_size

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["max_target_seqs"] == 100

    kwargs["hitlist_size"] = 50
    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["max_target_seqs"] == 50

    bad_int = [0,False,[],-1,None,str,"",{},int]
    for b in bad_int:
        print("passing bad hitlist_size:",b)
        kwargs["hitlist_size"] = b
        with pytest.raises(ValueError):
            sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)

    # -------------------------------------------------------------------------
    # evalue_cutoff

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["threshold"] == 0.01

    kwargs["e_value_cutoff"] = 0.001
    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["threshold"] == 0.001

    bad_float = [-1,[],None,str,"",{},int]
    for b in bad_float:

        kwargs["e_value_cutoff"] = b
        print("passing bad e_value_cutoff:",b)
        with pytest.raises(ValueError):
            sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)

    # -------------------------------------------------------------------------
    # gapcosts

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["gapopen"] == 11
    assert blast_kwargs["gapextend"] == 1

    kwargs["gapcosts"] = [10,1]
    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["gapopen"] == 10
    assert blast_kwargs["gapextend"] == 1

    bad_gap = [-1,[],[1,2,3],[-1,1],[1,-1],["test","this"],None,str,"",{},int]
    for b in bad_gap:
        kwargs["gapcosts"] = b
        print("passing bad gapcosts:",b)
        with pytest.raises(ValueError):
            sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)


    # -------------------------------------------------------------------------
    # Extra kwargs
    kwargs = copy.deepcopy(default_kwargs)
    kwargs["kwargs"] = {"extra_kwarg":7}
    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["extra_kwarg"] == 7


def test__construct_args(test_dataframes,tmpdir):

    # Make a fake blast db so code passes "file exists" check
    f = open(os.path.join(tmpdir,"GRCh38.psq"),"w")
    f.write("\n")
    f.close()
    fake_blast_db = os.path.join(tmpdir,"GRCh38")

    # Create some normal looking inputs to feed into _construct_args
    pfb_kwargs = get_public_param_defaults(local_blast,_pfb)
    pfb_kwargs["db"] = fake_blast_db
    pfb_kwargs["test_skip_blast_program_check"] = True
    pfb_kwargs["kwargs"] = {}
    df = test_dataframes["good-df"].copy()
    pfb_kwargs["sequence"] = df.sequence

    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**pfb_kwargs)

    # Run in configuration where we will have one query per core (5 threads,
    # 5 cores, 5 input sequences, long max_query_length)
    all_args, num_threads = _ca(sequence_list=sequence_list,
                                blast_function=blast_function,
                                blast_kwargs=blast_kwargs,
                                num_threads=5,
                                keep_tmp=False,
                                block_size=1,
                                test_num_cores=5)

    assert type(all_args) is list
    assert len(all_args) == 5
    assert num_threads == 5

    for i, a in enumerate(all_args):

        # Make sure it's pulling out sequences
        sequences = a[0]
        for j in range(len(sequences)):
            assert sequences[j] == df.sequence.iloc[j]

        # Make sure counter is working
        assert a[1][0] == i
        assert a[1][1] == i + 1

        assert a[2] == apps.NcbiblastpCommandline

        # useful kwargs
        assert a[3]["max_target_seqs"] == 100
        assert a[3]["threshold"] == 0.01
        assert a[3]["gapopen"] == 11
        assert a[3]["gapextend"] == 1

        # Num tries allowed
        assert a[4] == False

    # -------------------------------------------------------------------------
    # test sequence bits

    sequence_list = ["test"]
    all_args, num_threads = _ca(sequence_list=sequence_list,
                                blast_function=blast_function,
                                blast_kwargs=blast_kwargs,
                                num_threads=5,
                                keep_tmp=False,
                                block_size=1,
                                test_num_cores=5)

    assert len(all_args) == 1
    assert all_args[0][0][all_args[0][1][0]] == "test"
    assert num_threads == 1

    # Machine as two core, auto detect cores. Should have two args
    sequence_list = ["test","this"]
    all_args, num_threads = _ca(sequence_list=sequence_list,
                                blast_function=blast_function,
                                blast_kwargs=blast_kwargs,
                                num_threads=5,
                                keep_tmp=False,
                                block_size=1,
                                test_num_cores=2)

    assert len(all_args) == 2
    assert all_args[0][0][all_args[0][1][0]] == "test"
    assert all_args[1][0][all_args[1][1][0]] == "this"
    assert num_threads == 2

    # Machine as one core, auto detect cores. Should have one arg
    sequence_list = ["test","this"]
    all_args, num_threads = _ca(sequence_list,
                                blast_function=blast_function,
                                blast_kwargs=blast_kwargs,
                                block_size=20,
                                keep_tmp=False,
                                num_threads=-1,
                                test_num_cores=1)

    assert len(all_args) == 1
    assert all_args[0][0][all_args[0][1][0]] == "test"
    assert all_args[0][0][all_args[0][1][0]+1] == "this"
    assert num_threads == 1

    # -------------------------------------------------------------------------
    # block_size
    # making sure sequence bits are processed correctly when it's included

    sequence_list, blast_function, blast_kwargs, return_singleton = _pfb(**pfb_kwargs)

    bad_int = [0,False,[],-1,None,str,"",{},int]
    for b in bad_int:
        print("passing bad block_size:",b)
        with pytest.raises(ValueError):
            all_args, num_threads = _ca(sequence_list,
                                        blast_function=blast_function,
                                        blast_kwargs=blast_kwargs,
                                        block_size=b,
                                        keep_tmp=False,
                                        num_threads=-1,
                                        test_num_cores=1)




    all_args, num_threads = _ca(sequence_list,
                                blast_function=blast_function,
                                blast_kwargs=blast_kwargs,
                                block_size=5,
                                keep_tmp=False,
                                num_threads=-1,
                                test_num_cores=1)

    assert len(all_args) == 1
    assert all_args[0][1][0] == 0
    assert all_args[0][1][1] == 5

    # Make sure splitting looks reasonable -- each sequence on own
    all_args, num_threads = _ca(sequence_list,
                                blast_function=blast_function,
                                blast_kwargs=blast_kwargs,
                                block_size=1,
                                keep_tmp=False,
                                num_threads=-1,
                                test_num_cores=1)

    assert len(all_args) == 5
    for i, a in enumerate(all_args):
        seq = all_args[i][0][all_args[i][1][0]]
        assert seq == df.sequence.iloc[i]


    # Make sure splitting looks reasonable -- 2, 2, 1
    all_args, num_threads = _ca(sequence_list,
                                blast_function=blast_function,
                                blast_kwargs=blast_kwargs,
                                block_size=2,
                                keep_tmp=False,
                                num_threads=-1,
                                test_num_cores=1)

    assert len(all_args) == 3
    counter = 0
    for i, a in enumerate(all_args):
        seq = all_args[i][0][all_args[i][1][0]]
        assert seq == df.sequence.iloc[counter]
        counter += 2

    # Make sure splitting looks reasonable -- 3, 2
    all_args, num_threads = _ca(sequence_list,
                                blast_function=blast_function,
                                blast_kwargs=blast_kwargs,
                                block_size=3,
                                keep_tmp=False,
                                num_threads=-1,
                                test_num_cores=1)

    assert len(all_args) == 2
    counter = 0
    for i, a in enumerate(all_args):
        seq = all_args[i][0][all_args[i][1][0]]
        assert seq == df.sequence.iloc[counter]
        counter += 3

    # Make sure splitting looks reasonable -- 1
    all_args, num_threads = _ca(sequence_list,
                                blast_function=blast_function,
                                blast_kwargs=blast_kwargs,
                                block_size=5,
                                keep_tmp=False,
                                num_threads=-1,
                                test_num_cores=1)

    assert len(all_args) == 1

    # -------------------------------------------------------------------------
    # keep_tmp

    all_args, num_threads = _ca(sequence_list,
                                blast_function=blast_function,
                                blast_kwargs=blast_kwargs,
                                block_size=5,
                                keep_tmp=False,
                                num_threads=3,
                                test_num_cores=None)
    assert all_args[0][4] is False

    all_args, num_threads = _ca(sequence_list,
                                blast_function=blast_function,
                                blast_kwargs=blast_kwargs,
                                block_size=5,
                                keep_tmp=True,
                                num_threads=3,
                                test_num_cores=None)
    assert all_args[0][4] is True


    bad_bool = [1.5,[],None,str,"",{}]
    for b in bad_bool:
        print("passing bad keep_tmp:",b)
        with pytest.raises(ValueError):
            all_args, num_threads = _ca(sequence_list,
                                        blast_function=blast_function,
                                        blast_kwargs=blast_kwargs,
                                        block_size=5,
                                        keep_tmp=b,
                                        num_threads=3,
                                        test_num_cores=None)


    # -------------------------------------------------------------------------
    # num_threads.

    all_args, num_threads = _ca(sequence_list,
                                blast_function=blast_function,
                                blast_kwargs=blast_kwargs,
                                block_size=1,
                                keep_tmp=False,
                                num_threads=3,
                                test_num_cores=None)
    assert num_threads == 3


    bad_int = [0,False,[],-2,None,str,"",{},int]
    for b in bad_int:
        print("passing bad num_threads:",b)
        with pytest.raises(ValueError):
                    all_args, num_threads = _ca(sequence_list,
                                                blast_function=blast_function,
                                                blast_kwargs=blast_kwargs,
                                                block_size=5,
                                                keep_tmp=False,
                                                num_threads=b,
                                                test_num_cores=None)

def test__combine_hits(local_blast_output):

    # This is a set of two dataframes that have three and two outputs,
    # respectively. Make sure this is stitched into a reasonable set of outputs
    hits = copy.deepcopy(local_blast_output)

    df_list = _combine_hits(hits,return_singleton=False)
    assert type(df_list) is list
    assert len(df_list) == 5
    for i in range(len(df_list)):
        assert df_list[i]["query"].iloc[i] == f"count{i}"

    single_mask = hits[0]["query"] == "count0"
    single_hit = [hits[0].loc[single_mask,:]]

    df_list = _combine_hits(single_hit,return_singleton=True)
    assert type(df_list) is type(pd.DataFrame({"test":[1,2,3]}))

    df_list = _combine_hits(single_hit,return_singleton=False)
    assert type(df_list) is list
    assert len(df_list) == 1
