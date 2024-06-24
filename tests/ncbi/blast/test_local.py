
import pytest
from conftest import get_public_param_defaults

import topiary
from topiary.ncbi.blast.local import local_blast
from topiary.ncbi.blast.local import _prepare_for_blast
from topiary.ncbi.blast.local import _construct_args
from topiary.ncbi.blast.local import _combine_hits
from topiary.ncbi.blast.local import _local_blast_thread_function

# Functions used in testing
from topiary.ncbi.blast.make import make_blast_db

import numpy as np
import pandas as pd

import copy
import os
import glob

def _blast_args_to_keys(blast_args):
    """
    Take a list of args like ["blastp","-db","stupid",...] and return a 
    dictionary where every arg is a key in a dictionary pointing to a value.
    Assumes a single value for each key and ignores anything not starting 
    with -. Would return {'-db':'stupid'} from above.
    """
    
    out_dict = {}
    i = 0
    while i < len(blast_args):
        if blast_args[i].startswith('-'):
            out_dict[blast_args[i]] = blast_args[i+1]
            i += 1
        i += 1

    return out_dict      

def test__prepare_for_blast(test_dataframes,tmpdir):

    # Make a fake blast db so code passes "file exists" check
    f = open(os.path.join(tmpdir,"GRCh38.psq"),"w")
    f.write("\n")
    f.close()
    fake_blast_db = os.path.join(tmpdir,"GRCh38")

    default_kwargs = get_public_param_defaults(local_blast,_prepare_for_blast)
    default_kwargs["db"] = fake_blast_db
    default_kwargs["test_skip_blast_program_check"] = True
    default_kwargs["kwargs"] = {}
    df = test_dataframes["good-df"].copy()
    default_kwargs["sequence"] = df.sequence

    # Skip the check for a working blast program.
    sequence_list, blast_args, return_singleton = _prepare_for_blast(**default_kwargs)

    # Should be a list of length 1 holding a single dictionary
    assert type(sequence_list) is list
    assert len(sequence_list) == 5
    assert type(blast_args) is list
    assert return_singleton == False

    # Make sure it got the sequences right
    assert np.array_equal(sequence_list,df.sequence)

    # Make sure the blast command is constructed correctly
    assert blast_args[0] == "blastp"
    assert blast_args[1] == "-db"
    assert os.path.split(blast_args[2])[-1] == "GRCh38"
    
    expected = ['-outfmt', '5',
                '-max_target_seqs','100',
                '-threshold', '1.00000e-03',
                '-gapopen', '11',
                '-gapextend', '1']
    for i in range(len(expected)):
        assert blast_args[i+3] == expected[i]

    # -------------------------------------------------------------------------
    # test sequence bits

    kwargs = copy.deepcopy(default_kwargs)

    kwargs["sequence"] = "test"
    sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)
    assert len(sequence_list) == 1
    assert sequence_list[0] == "test"
    assert return_singleton == True

    kwargs["sequence"] = ["test"]
    sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)
    assert len(sequence_list) == 1
    assert sequence_list[0] == "test"
    assert return_singleton == False

    kwargs["sequence"] = ["test","this"]
    sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)
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
            sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)

        kwargs["sequence"] = [b]
        with pytest.raises(ValueError):
            sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)

    # -------------------------------------------------------------------------
    # test db

    kwargs = copy.deepcopy(default_kwargs)

    kwargs["db"] = "not_really_db"
    with pytest.raises(FileNotFoundError):
        sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)

    bad_db = [1,False,1.0,-1,None,str,"",{}]
    for b in bad_db:
        print("passing bad db:",b)
        kwargs["db"] = b
        with pytest.raises(FileNotFoundError):
            sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)

    # -------------------------------------------------------------------------
    # blast_program

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)
    assert blast_args[0] == "blastp"

    kwargs["blast_program"] = "tblastx"
    sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)
    assert  blast_args[0] == "tblastx"

    bad_pg = ["not_really_program",1,False,1.0,-1,None,str,"",{}]
    for b in bad_pg:
        print("passing bad blast_program:",b)
        kwargs["blast_program"] = b
        with pytest.raises(ValueError):
            sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)


    # -------------------------------------------------------------------------
    # hitlist_size

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)
    blast_kwargs = _blast_args_to_keys(blast_args)

    assert blast_kwargs["-max_target_seqs"] == "100"

    kwargs["hitlist_size"] = "50"
    sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)
    blast_kwargs = _blast_args_to_keys(blast_args)
    assert blast_kwargs["-max_target_seqs"] == "50"

    bad_int = [0,False,[],-1,None,str,"",{},int]
    for b in bad_int:
        print("passing bad hitlist_size:",b)
        kwargs["hitlist_size"] = b
        with pytest.raises(ValueError):
            sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)


    # -------------------------------------------------------------------------
    # evalue_cutoff

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)
    assert blast_kwargs["-threshold"] == "1.00000e-03"

    kwargs["e_value_cutoff"] = 0.01
    sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)
    blast_kwargs = _blast_args_to_keys(blast_args)
    assert blast_kwargs["-threshold"] == "1.00000e-02" 

    bad_float = [-1,[],None,str,"",{},int]
    for b in bad_float:

        kwargs["e_value_cutoff"] = b
        print("passing bad e_value_cutoff:",b)
        with pytest.raises(ValueError):
            sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)

    # -------------------------------------------------------------------------
    # gapcosts

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)
    blast_kwargs = _blast_args_to_keys(blast_args)
    assert blast_kwargs["-gapopen"] == "11"
    assert blast_kwargs["-gapextend"] == "1"

    kwargs["gapcosts"] = [10,1]
    sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)
    blast_kwargs = _blast_args_to_keys(blast_args)
    assert blast_kwargs["-gapopen"] == "10"
    assert blast_kwargs["-gapextend"] == "1"

    bad_gap = [-1,[],[1,2,3],[-1,1],[1,-1],["test","this"],None,str,"",{},int]
    for b in bad_gap:
        kwargs["gapcosts"] = b
        print("passing bad gapcosts:",b)
        with pytest.raises(ValueError):
            sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)

    # -------------------------------------------------------------------------
    # Extra kwargs
    kwargs = copy.deepcopy(default_kwargs)
    kwargs["kwargs"] = {"extra_kwarg":7}
    sequence_list, blast_args, return_singleton = _prepare_for_blast(**kwargs)
    blast_kwargs = _blast_args_to_keys(blast_args)
    assert blast_kwargs["-extra_kwarg"] == "7"

def test__construct_args(test_dataframes,tmpdir):

    # Make a fake blast db so code passes "file exists" check
    f = open(os.path.join(tmpdir,"GRCh38.psq"),"w")
    f.write("\n")
    f.close()
    fake_blast_db = os.path.join(tmpdir,"GRCh38")

    # Create some normal looking inputs to feed into _construct_args
    pfb_kwargs = get_public_param_defaults(local_blast,_prepare_for_blast)
    pfb_kwargs["db"] = fake_blast_db
    pfb_kwargs["test_skip_blast_program_check"] = True
    pfb_kwargs["kwargs"] = {}
    df = test_dataframes["good-df"].copy()
    pfb_kwargs["sequence"] = df.sequence

    sequence_list, blast_args, return_singleton = _prepare_for_blast(**pfb_kwargs)

    # Run in configuration where we will have one query per core (5 threads,
    # 5 cores, 5 input sequences, long max_query_length)
    kwargs_list, num_threads = _construct_args(sequence_list=sequence_list,
                                   blast_args=blast_args,
                                   num_threads=5,
                                   keep_blast_xml=False,
                                   block_size=1,
                                   manual_num_cores=5)

    assert type(kwargs_list) is list
    assert len(kwargs_list) == 5
    assert num_threads == 5

    for i, a in enumerate(kwargs_list):

        # Make sure it's pulling out sequences
        sequences = a["sequence_list"]
        for j in range(len(sequences)):
            assert sequences[j] == df.sequence.iloc[j]

        # Make sure counter is working
        assert a["index"][0] == i
        assert a["index"][1] == i + 1
        

        blast_kwargs = _blast_args_to_keys(a["blast_args"])

        # useful kwargs
        assert a["blast_args"][0] == "blastp"
        assert blast_kwargs["-max_target_seqs"] == "100"
        assert blast_kwargs["-threshold"] == "1.00000e-03"
        assert blast_kwargs["-gapopen"] == "11"
        assert blast_kwargs["-gapextend"] == "1"

        # Keep tmp
        assert a["keep_blast_xml"] == False

    # -------------------------------------------------------------------------
    # test sequence bits

    sequence_list = ["test"]
    all_args, num_threads = _construct_args(sequence_list=sequence_list,
                                blast_args=blast_args,
                                num_threads=5,
                                keep_blast_xml=False,
                                block_size=1,
                                manual_num_cores=5)

    assert len(all_args) == 1
    assert all_args[0]["sequence_list"][all_args[0]["index"][0]] == "test"
    assert num_threads == 1

    # Machine as two core, auto detect cores. Should have two args
    sequence_list = ["test","this"]
    all_args, num_threads = _construct_args(sequence_list=sequence_list,
                                blast_args=blast_args,
                                num_threads=5,
                                keep_blast_xml=False,
                                block_size=1,
                                manual_num_cores=2)

    assert len(all_args) == 2
    assert all_args[0]["sequence_list"][all_args[0]["index"][0]] == "test"
    assert all_args[1]["sequence_list"][all_args[1]["index"][0]] == "this"
    assert num_threads == 2

    # Machine as one core, auto detect cores. Should have one arg
    sequence_list = ["test","this"]
    all_args, num_threads = _construct_args(sequence_list,
                                blast_args=blast_args,
                                block_size=20,
                                keep_blast_xml=False,
                                num_threads=-1,
                                manual_num_cores=1)

    assert len(all_args) == 1
    assert all_args[0]["sequence_list"][all_args[0]["index"][0]] == "test"
    assert all_args[0]["sequence_list"][all_args[0]["index"][0]+1] == "this"
    assert num_threads == 1

    # -------------------------------------------------------------------------
    # block_size
    # making sure sequence bits are processed correctly when it's included

    sequence_list, blast_args, return_singleton = _prepare_for_blast(**pfb_kwargs)

    bad_int = [0,False,[],-1,None,str,"",{},int]
    for b in bad_int:
        print("passing bad block_size:",b)
        with pytest.raises(ValueError):
            all_args, num_threads = _construct_args(sequence_list,
                                        blast_args=blast_args,
                                        block_size=b,
                                        keep_blast_xml=False,
                                        num_threads=-1,
                                        manual_num_cores=1)




    all_args, num_threads = _construct_args(sequence_list,
                                blast_args=blast_args,
                                block_size=5,
                                keep_blast_xml=False,
                                num_threads=-1,
                                manual_num_cores=1)

    assert len(all_args) == 1
    assert all_args[0]["index"][0] == 0
    assert all_args[0]["index"][1] == 5

    # Make sure splitting looks reasonable -- each sequence on own
    all_args, num_threads = _construct_args(sequence_list,
                                blast_args=blast_args,
                                block_size=1,
                                keep_blast_xml=False,
                                num_threads=-1,
                                manual_num_cores=1)

    assert len(all_args) == 5
    for i, a in enumerate(all_args):
        seq = all_args[i]["sequence_list"][all_args[i]["index"][0]]
        assert seq == df.sequence.iloc[i]


    # Make sure splitting looks reasonable -- 2, 2, 1
    all_args, num_threads = _construct_args(sequence_list,
                                blast_args=blast_args,
                                block_size=2,
                                keep_blast_xml=False,
                                num_threads=-1,
                                manual_num_cores=1)

    assert len(all_args) == 3
    counter = 0
    for i, a in enumerate(all_args):
        seq = all_args[i]["sequence_list"][all_args[i]["index"][0]]
        assert seq == df.sequence.iloc[counter]
        counter += 2

    # Make sure splitting looks reasonable -- 3, 2
    all_args, num_threads = _construct_args(sequence_list,
                                blast_args=blast_args,
                                block_size=3,
                                keep_blast_xml=False,
                                num_threads=-1,
                                manual_num_cores=1)

    assert len(all_args) == 2
    counter = 0
    for i, a in enumerate(all_args):
        seq = all_args[i]["sequence_list"][all_args[i]["index"][0]]
        assert seq == df.sequence.iloc[counter]
        counter += 3

    # Make sure splitting looks reasonable -- 1
    all_args, num_threads = _construct_args(sequence_list,
                                blast_args=blast_args,
                                block_size=5,
                                keep_blast_xml=False,
                                num_threads=-1,
                                manual_num_cores=1)

    assert len(all_args) == 1

    # -------------------------------------------------------------------------
    # keep_blast_xml

    all_args, num_threads = _construct_args(sequence_list,
                                blast_args=blast_args,
                                block_size=5,
                                keep_blast_xml=False,
                                num_threads=3,
                                manual_num_cores=None)
    assert all_args[0]["keep_blast_xml"] is False

    all_args, num_threads = _construct_args(sequence_list,
                                blast_args=blast_args,
                                block_size=5,
                                keep_blast_xml=True,
                                num_threads=3,
                                manual_num_cores=None)
    assert all_args[0]["keep_blast_xml"] is True


    bad_bool = [1.5,[],None,str,"",{}]
    for b in bad_bool:
        print("passing bad keep_blast_xml:",b)
        with pytest.raises(ValueError):
            all_args, num_threads = _construct_args(sequence_list,
                                        blast_args=blast_args,
                                        block_size=5,
                                        keep_blast_xml=b,
                                        num_threads=3,
                                        manual_num_cores=None)


    # -------------------------------------------------------------------------
    # num_threads.

    all_args, num_threads = _construct_args(sequence_list,
                                blast_args=blast_args,
                                block_size=1,
                                keep_blast_xml=False,
                                num_threads=3,
                                manual_num_cores=1000)
    assert num_threads == 3


    bad_int = [0,False,[],-2,None,str,"",{},int]
    for b in bad_int:
        print("passing bad num_threads:",b)
        with pytest.raises(ValueError):
            all_args, num_threads = _construct_args(sequence_list,
                                        blast_args=blast_args,
                                        block_size=5,
                                        keep_blast_xml=False,
                                        num_threads=b,
                                        manual_num_cores=None)

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

@pytest.mark.skipif(os.name == "nt",reason="blast cannot be installed via conda on windows")
@pytest.mark.run_blast
def test__local_blast_thread_function(test_dataframes,
                                      make_blast_db_files,
                                      tmpdir):

    # Move into temporary directory
    cwd = os.getcwd()
    os.chdir(tmpdir)

    # Make a blast database
    faa = [make_blast_db_files["test1.faa"],make_blast_db_files["test2.faa"]]
    out = "blastdb"
    make_blast_db(faa,db_name=out)

    pfb_kwargs = get_public_param_defaults(local_blast,_prepare_for_blast)
    pfb_kwargs["db"] = "blastdb"
    pfb_kwargs["test_skip_blast_program_check"] = False
    pfb_kwargs["kwargs"] = {}
    df = test_dataframes["good-df"].copy()
    pfb_kwargs["sequence"] = df.sequence

    sequence_list, blast_args, _ = _prepare_for_blast(**pfb_kwargs)

    out_df = _local_blast_thread_function(sequence_list=sequence_list,
                                          index=[0,len(sequence_list)],
                                          blast_args=blast_args,
                                          keep_blast_xml=True)

    # Should be five sequences and should have created an xml file
    assert len(out_df) == 5
    xml_files = glob.glob("*blast-out.xml") 
    assert len(xml_files) == 1
    for f in xml_files:
        os.remove(f)

    # Should be two sequences, single xml file
    out_df = _local_blast_thread_function(sequence_list=sequence_list,
                                          index=[0,2],
                                          blast_args=blast_args,
                                          keep_blast_xml=True)
    
    assert len(out_df) == 2
    xml_files = glob.glob("*blast-out.xml") 
    assert len(xml_files) == 1
    for f in xml_files:
        os.remove(f)

    os.chdir(cwd)

@pytest.mark.skipif(os.name == "nt",reason="blast cannot be installed via conda on windows")
@pytest.mark.run_blast
def test_local_blast(test_dataframes,
                     make_blast_db_files,
                     tmpdir):

    # This is a fairly minimal test suite because this function simply wraps and
    # the set of local functions that are tested extensively in the rest of this
    # test file. 

    # Move into temporary directory
    cwd = os.getcwd()
    os.chdir(tmpdir)

    # Make a blast database
    faa = [make_blast_db_files["test1.faa"],make_blast_db_files["test2.faa"]]
    out = "blastdb"
    make_blast_db(faa,db_name=out)

    df = test_dataframes["good-df"].copy()
    sequence = df.sequence

    query_hits = local_blast(sequence=sequence,
                             db="blastdb")

    assert len(sequence) == len(query_hits)

    os.chdir(cwd)
