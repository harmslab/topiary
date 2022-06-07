
import pytest
from conftest import get_public_param_defaults

import topiary
from topiary.external.ncbi._ncbi_blast import ncbi_blast as ncbi_blast
from topiary.external.ncbi._ncbi_blast import _prepare_for_blast as _pfb
from topiary.external.ncbi._ncbi_blast import _construct_args as _ca
from topiary.external.ncbi._ncbi_blast import _combine_hits

import numpy as np
import pandas as pd

import copy
import multiprocessing as mp

def test__prepare_blast(test_dataframes):

    default_kwargs = get_public_param_defaults(ncbi_blast,_pfb)
    default_kwargs["kwargs"] = {}
    df = test_dataframes["good-df"].copy()
    default_kwargs = copy.deepcopy(default_kwargs)
    default_kwargs["sequence"] = df.sequence

    sequence_list, blast_kwargs, return_singleton = _pfb(**default_kwargs)

    # Should be a list of length 1 holding a single dictionary
    assert type(sequence_list) is list
    assert len(sequence_list) == 5
    assert type(blast_kwargs) is dict
    assert return_singleton == False

    # Make sure it got the sequences right
    assert np.array_equal(sequence_list,df.sequence)

    # Make sure the rest of the arguments were done properly
    assert blast_kwargs["database"] == "nr"
    assert blast_kwargs["hitlist_size"] == '50'
    assert blast_kwargs["program"] == "blastp"
    assert blast_kwargs["expect"] == '0.01'
    assert blast_kwargs["gapcosts"] == '11 1'
    assert blast_kwargs["url_base"] == "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    with pytest.raises(KeyError):
        blast_kwargs["entrez_query"]

    # -------------------------------------------------------------------------
    # test sequence bits

    kwargs = copy.deepcopy(default_kwargs)

    kwargs["sequence"] = "test"
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert len(sequence_list) == 1
    assert sequence_list[0] == "test"
    assert return_singleton == True

    kwargs["sequence"] = ["test"]
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert len(sequence_list) == 1
    assert sequence_list[0] == "test"
    assert return_singleton == False

    kwargs["sequence"] = ["test","this"]
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
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
            sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)

        kwargs["sequence"] = [b]
        with pytest.raises(ValueError):
            sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)

    # -------------------------------------------------------------------------
    # test db

    kwargs = copy.deepcopy(default_kwargs)

    kwargs["sequence"] = df.sequence
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["database"] == "nr"

    kwargs["db"] = "not_really_db"
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["database"] == "not_really_db"

    bad_db = [1,False,1.0,-1,None,str,"",{}]
    for b in bad_db:
        print("passing bad db:",b)
        kwargs["db"] = b
        with pytest.raises(ValueError):
            sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)

    # -------------------------------------------------------------------------
    # taxid

    kwargs = copy.deepcopy(default_kwargs)

    kwargs["taxid"] = 9606
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["entrez_query"] == "txid9606[ORGN]"

    kwargs["taxid"] = "9606"
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["entrez_query"] == "txid9606[ORGN]"

    kwargs["taxid"] = (9606,1234)
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["entrez_query"] == "txid9606[ORGN] or txid1234[ORGN]"

    kwargs["taxid"] = ["9606","1234"]
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["entrez_query"] == "txid9606[ORGN] or txid1234[ORGN]"

    kwargs["taxid"] = None
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    with pytest.raises(KeyError):
        blast_kwargs["entrez_query"]

    bad_taxid = [1.15,str,int,[None,None]]
    for b in bad_taxid:

        kwargs["taxid"] = b
        with pytest.raises(ValueError):
            print("passing bad taxid:",b)
            sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
        kwargs["taxid"] = [9606,b]
        with pytest.raises(ValueError):
            print("passing bad taxid:",[9606,b])
            sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)

    # -------------------------------------------------------------------------
    # blast_program

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["program"] == "blastp"

    kwargs["blast_program"] = "not_really_pg"
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["program"] == "not_really_pg"

    bad_pg = [1,False,1.0,-1,None,str,"",{}]
    for b in bad_pg:
        print("passing bad blast_program:",b)
        kwargs["blast_program"] = b
        with pytest.raises(ValueError):
            sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)

    # -------------------------------------------------------------------------
    # histlist_size

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["hitlist_size"] == '50'

    kwargs["hitlist_size"] = 100
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["hitlist_size"] == '100'

    bad_int = [0,False,[],-1,None,str,"",{},int]
    for b in bad_int:
        print("passing bad hitlist_size:",b)
        kwargs["hitlist_size"] = b
        with pytest.raises(ValueError):
            sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)



    # -------------------------------------------------------------------------
    # evalue_cutoff

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["expect"] == '0.01'

    kwargs["e_value_cutoff"] = 0.001
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["expect"] == '0.001'

    bad_float = [-1,[],None,str,"",{},int]
    for b in bad_float:

        kwargs["e_value_cutoff"] = b
        print("passing bad e_value_cutoff:",b)
        with pytest.raises(ValueError):
            sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)

    # -------------------------------------------------------------------------
    # gapcosts

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["gapcosts"] == '11 1'

    kwargs["gapcosts"] = [10,1]
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["gapcosts"] == '10 1'

    bad_gap = [-1,[],[1,2,3],[-1,1],[1,-1],["test","this"],None,str,"",{},int]
    for b in bad_gap:
        kwargs["gapcosts"] = b
        print("passing bad gapcosts:",b)
        with pytest.raises(ValueError):
            sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)

    # -------------------------------------------------------------------------
    # url_base
    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["url_base"] == "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

    kwargs["url_base"] = "not_really_pg"
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["url_base"] == "not_really_pg"

    bad_pg = [1,False,1.0,-1,None,str,"",{}]
    for b in bad_pg:
        print("passing bad url_base:",b)
        kwargs["url_base"] = b
        with pytest.raises(ValueError):
            sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)

    # -------------------------------------------------------------------------
    # Extra kwargs
    kwargs = copy.deepcopy(default_kwargs)
    kwargs["kwargs"] = {"extra_kwarg":7}
    sequence_list, blast_kwargs, return_singleton = _pfb(**kwargs)
    assert blast_kwargs["extra_kwarg"] == 7

def test__construct_args(test_dataframes):

    # Create some normal looking inputs to feed into _construct_args
    pfb_kwargs = get_public_param_defaults(ncbi_blast,_pfb)
    pfb_kwargs["kwargs"] = {}
    df = test_dataframes["good-df"].copy()
    pfb_kwargs["sequence"] = df.sequence

    sequence_list, blast_kwargs, return_singleton = _pfb(**pfb_kwargs)

    # Run in configuration where we will have one query per core (5 threads,
    # 5 cores, 5 input sequences, long max_query_length)
    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=10000,
                                num_tries_allowed=5,
                                num_threads=5,
                                test_num_cores=5)

    assert type(all_args) is list
    assert len(all_args) == 5
    assert num_threads == 5

    for i, a in enumerate(all_args):

        # Make sure it's pulling out sequences
        seq = a[0]["sequence"].split("\n")[1].strip()
        assert seq == df.sequence.iloc[i]

        # useful kwargs
        assert a[0]["database"] == "nr"
        assert a[0]["hitlist_size"] == '50'
        assert a[0]["program"] == "blastp"
        assert a[0]["expect"] == '0.01'
        assert a[0]["gapcosts"] == '11 1'
        assert a[0]["url_base"] == "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

        # Make sure counter is working
        assert a[1] == i

        # Num tries allowed
        assert a[2] == 5

    # -------------------------------------------------------------------------
    # test sequence bits

    sequence_list = ["test"]
    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=10000,
                                num_tries_allowed=5,
                                num_threads=5,
                                test_num_cores=5)

    assert len(all_args) == 1
    assert all_args[0][0]["sequence"].split("\n")[1] == "test"
    assert num_threads == 1

    # Machine as two core, auto detect cores. Should have two args
    sequence_list = ["test","this"]
    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=10000,
                                num_tries_allowed=5,
                                num_threads=-1,
                                test_num_cores=2)

    assert len(all_args) == 2
    assert all_args[0][0]["sequence"].split("\n")[1] == "test"
    assert all_args[1][0]["sequence"].split("\n")[1] == "this"
    assert num_threads == 2

    # Machine as one core, auto detect cores. Should have one arg
    sequence_list = ["test","this"]
    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=10000,
                                num_tries_allowed=5,
                                num_threads=-1,
                                test_num_cores=1)

    assert len(all_args) == 1
    assert all_args[0][0]["sequence"].split("\n")[1] == "test"
    assert all_args[0][0]["sequence"].split("\n")[3] == "this"
    assert num_threads == 1

    # Machine as one core. Pass in 2. Should have one arg, one thread
    sequence_list = ["test","this"]
    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=10000,
                                num_tries_allowed=5,
                                num_threads=2,
                                test_num_cores=1)

    assert len(all_args) == 1
    assert all_args[0][0]["sequence"].split("\n")[1] == "test"
    assert all_args[0][0]["sequence"].split("\n")[3] == "this"
    assert num_threads == 1


    # One core, but super short query length. Should give two args, but one
    # thread
    sequence_list = [25*"test",25*"this"]
    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=150,
                                num_tries_allowed=5,
                                num_threads=-1,
                                test_num_cores=1)
    assert len(all_args) == 2
    assert all_args[0][0]["sequence"].split("\n")[1] == 25*"test"
    assert all_args[1][0]["sequence"].split("\n")[1] == 25*"this"
    assert num_threads == 1

    # -------------------------------------------------------------------------
    # max_query_length
    # making sure sequence bits are processed correctly when it's included

    sequence_list, blast_kwargs, return_singleton = _pfb(**pfb_kwargs)

    bad_int = [0,False,[],-1,None,str,"",{},int]
    for b in bad_int:
        print("passing bad max_query_length:",b)
        with pytest.raises(ValueError):
            all_args, num_threads = _ca(sequence_list,
                                        blast_kwargs=blast_kwargs,
                                        max_query_length=b,
                                        num_tries_allowed=5,
                                        num_threads=-1,
                                        test_num_cores=1)

    # Get expected and actual sequence length for this df.sequence compiled
    # into >countX\nSEQUENCE\n ... format. Assumes there are less than 10
    # seqs in test dataset. Also assumes that the test sequences all have the
    # same length (160 amino acids)
    lens = [len(s) for s in df.sequence]
    expected_length = np.sum(lens) + len(lens)*len(">countX\n\n") - 1

    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=10000,
                                num_tries_allowed=5,
                                num_threads=1,
                                test_num_cores=1)

    assert len(all_args) == 1
    seqs = all_args[0][0]["sequence"].split("\n")
    idx = [1,3,5,7,9]
    for i in range(5):
        assert seqs[idx[i]].strip() == df.sequence.iloc[i]

    # If our max query length is shorter than one of our sequences, it should
    # catch and die
    long_indiv_sequence = np.max([len(s) for s in df.sequence])
    with pytest.raises(ValueError):
        all_args, num_threads = _ca(sequence_list,
                                    blast_kwargs=blast_kwargs,
                                    max_query_length=long_indiv_sequence//2,
                                    num_tries_allowed=5,
                                    num_threads=1,
                                    test_num_cores=1)

    # Make sure splitting looks reasonable -- each sequence on own
    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=180,
                                num_tries_allowed=5,
                                num_threads=1,
                                test_num_cores=1)

    assert len(all_args) == 5
    for i, a in enumerate(all_args):

        # Make sure it's pulling out sequences
        seq = a[0]["sequence"].split("\n")[1].strip()
        assert seq == df.sequence.iloc[i]


    # Make sure splitting looks reasonable -- 2, 2, 1
    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=340,
                                num_tries_allowed=5,
                                num_threads=1,
                                test_num_cores=1)

    assert len(all_args) == 3
    counter = 0
    for i, a in enumerate(all_args):

        # Make sure it's pulling out sequences
        seq = a[0]["sequence"].split("\n")[1].strip()
        assert seq == df.sequence.iloc[counter]
        counter += 2

    # Make sure splitting looks reasonable -- 3, 2
    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=510,
                                num_tries_allowed=5,
                                num_threads=1,
                                test_num_cores=1)

    assert len(all_args) == 2
    counter = 0
    for i, a in enumerate(all_args):

        # Make sure it's pulling out sequences
        seq = a[0]["sequence"].split("\n")[1].strip()
        assert seq == df.sequence.iloc[counter]
        counter += 3


    # Make sure splitting looks reasonable -- 1
    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=10000,
                                num_tries_allowed=5,
                                num_threads=1,
                                test_num_cores=1)

    assert len(all_args) == 1

    # -------------------------------------------------------------------------
    # num_tries_allowed.

    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=10000,
                                num_tries_allowed=7,
                                num_threads=1,
                                test_num_cores=1)
    assert all_args[0][2] == 7


    bad_int = [0,False,[],-1,None,str,"",{},int]
    for b in bad_int:
        print("passing bad num_tries_allowed:",b)
        with pytest.raises(ValueError):
                all_args, num_threads = _ca(sequence_list,
                                            blast_kwargs=blast_kwargs,
                                            max_query_length=10000,
                                            num_tries_allowed=b,
                                            num_threads=1,
                                            test_num_cores=1)

    # -------------------------------------------------------------------------
    # num_threads.

    all_args, num_threads = _ca(sequence_list,
                                blast_kwargs=blast_kwargs,
                                max_query_length=10000,
                                num_tries_allowed=7,
                                num_threads=3,
                                test_num_cores=None)
    assert num_threads == 3


    bad_int = [0,False,[],-2,None,str,"",{},int]
    for b in bad_int:
        print("passing bad num_threads:",b)
        with pytest.raises(ValueError):
                all_args, num_threads = _ca(sequence_list,
                                            blast_kwargs=blast_kwargs,
                                            max_query_length=10000,
                                            num_tries_allowed=5,
                                            num_threads=b,
                                            test_num_cores=None)

def test__combine_hits(ncbi_blast_server_output):

    # This is a set of two dataframes that have three and two outputs,
    # respectively. Make sure this is stitched into a reasonable set of outputs
    hits = copy.deepcopy(ncbi_blast_server_output)

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
