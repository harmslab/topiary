
import pytest
from topiary.ncbi import ncbi_blast

import numpy as np

def test_ncbi_blast(test_dataframes):


    df = test_dataframes["good-df"]
    out = ncbi_blast(sequence=df.sequence,test_run=True)

    # Should be a list of length 1 holding a single dictionary
    assert type(out) is list
    assert len(out) == 1
    assert type(out[0]) is dict

    # Make sure it got the sequences right
    blast_kwargs = out[0]
    sequences = blast_kwargs["sequence"]
    sequences = sequences.split("\n")
    sequences = [sequences[i] for i in range(1,len(sequences),2)]

    assert np.array_equal(sequences,df.sequence)

    # Make sure other arguments are processing properly
    out = ncbi_blast(sequence=df.sequence,test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["database"] == "nr"
    assert blast_kwargs["hitlist_size"] == '50'
    assert blast_kwargs["program"] == "blastp"
    assert blast_kwargs["expect"] == '0.01'
    assert blast_kwargs["gapcosts"] == '11 1'
    assert blast_kwargs["url_base"] == "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    with pytest.raises(KeyError):
        blast_kwargs["entrez_query"]

    # test sequence bits
    out = ncbi_blast(sequence="test",test_run=True)
    blast_kwargs = out[0]
    sequences = blast_kwargs["sequence"]
    assert sequences.split("\n")[1] == "test"

    out = ncbi_blast(sequence=["test"],test_run=True)
    blast_kwargs = out[0]
    sequences = blast_kwargs["sequence"]
    assert sequences.split("\n")[1] == "test"

    out = ncbi_blast(sequence=["test","this"],test_run=True)
    blast_kwargs = out[0]
    sequences = blast_kwargs["sequence"]
    assert sequences.split("\n")[1] == "test"
    assert sequences.split("\n")[3] == "this"

    # Send in bad sequences
    bad_sequences = [1,False,1.0,-1,None,str,""]
    for b in bad_sequences:
        print("passing bad sequence:",b)
        with pytest.raises(ValueError):
            out = ncbi_blast(sequence=b,test_run=True)
        with pytest.raises(ValueError):
            out = ncbi_blast(sequence=[b],test_run=True)

    # db
    out = ncbi_blast(sequence=df.sequence,test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["database"] == "nr"

    out = ncbi_blast(sequence=df.sequence,db="not_really_db",test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["database"] == "not_really_db"

    bad_db = [1,False,1.0,-1,None,str,"",{}]
    for b in bad_db:
        print("passing bad db:",b)
        with pytest.raises(ValueError):
            out = ncbi_blast(sequence=df.sequence,db=b,test_run=True)

    # taxid
    out = ncbi_blast(sequence=df.sequence,taxid=9606,test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["entrez_query"] == "txid9606[ORGN]"

    out = ncbi_blast(sequence=df.sequence,taxid="9606",test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["entrez_query"] == "txid9606[ORGN]"

    out = ncbi_blast(sequence=df.sequence,taxid=(9606,1234),test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["entrez_query"] == "txid9606[ORGN] or txid1234[ORGN]"

    out = ncbi_blast(sequence=df.sequence,taxid=["9606","1234"],test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["entrez_query"] == "txid9606[ORGN] or txid1234[ORGN]"

    out = ncbi_blast(sequence=df.sequence,taxid=None,test_run=True)
    blast_kwargs = out[0]
    with pytest.raises(KeyError):
        blast_kwargs["entrez_query"]

    bad_taxid = [1.15,str,int,[None,None]]
    for b in bad_taxid:
        with pytest.raises(ValueError):
            print("passing bad taxid:",b)
            out = ncbi_blast(sequence=df.sequence,taxid=b,test_run=True)
        with pytest.raises(ValueError):
            print("passing bad taxid:",[9606,b])
            out = ncbi_blast(sequence=df.sequence,taxid=[9606,b],test_run=True)

    # blast_program
    out = ncbi_blast(sequence=df.sequence,test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["program"] == "blastp"

    out = ncbi_blast(sequence=df.sequence,blast_program="not_really_pg",test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["program"] == "not_really_pg"

    bad_pg = [1,False,1.0,-1,None,str,"",{}]
    for b in bad_pg:
        print("passing bad blast_program:",b)
        with pytest.raises(ValueError):
            out = ncbi_blast(sequence=df.sequence,blast_program=b,test_run=True)

    # hitlist_size
    out = ncbi_blast(sequence=df.sequence,test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["hitlist_size"] == '50'

    out = ncbi_blast(sequence=df.sequence,hitlist_size=100,test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["hitlist_size"] == '100'

    bad_int = [0,False,[],-1,None,str,"",{},int]
    for b in bad_int:
        print("passing bad hitlist_size:",b)
        with pytest.raises(ValueError):
            out = ncbi_blast(sequence=df.sequence,hitlist_size=b,test_run=True)

    # e_value_cutoff
    out = ncbi_blast(sequence=df.sequence,test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["expect"] == '0.01'

    out = ncbi_blast(sequence=df.sequence,e_value_cutoff=0.001,test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["expect"] == '0.001'

    bad_float = [-1,[],None,str,"",{},int]
    for b in bad_float:
        print("passing bad e_value_cutoff:",b)
        with pytest.raises(ValueError):
            out = ncbi_blast(sequence=df.sequence,e_value_cutoff=b,test_run=True)

    # gapcosts
    out = ncbi_blast(sequence=df.sequence,test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["gapcosts"] == '11 1'

    out = ncbi_blast(sequence=df.sequence,gapcosts=[10,1],test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["gapcosts"] == '10 1'

    bad_gap = [-1,[],[1,2,3],[-1,1],[1,-1],["test","this"],None,str,"",{},int]
    for b in bad_gap:
        print("passing bad gapcosts:",b)
        with pytest.raises(ValueError):
            out = ncbi_blast(sequence=df.sequence,gapcosts=b,test_run=True)

    # num_tries_allowed. can't check output, but can at least make sure error
    # checking is working
    out = ncbi_blast(sequence=df.sequence,test_run=True)
    out = ncbi_blast(sequence=df.sequence,num_tries_allowed=100,test_run=True)
    bad_int = [0,False,[],-1,None,str,"",{},int]
    for b in bad_int:
        print("passing bad num_tries_allowed:",b)
        with pytest.raises(ValueError):
            out = ncbi_blast(sequence=df.sequence,num_tries_allowed=b,test_run=True)

    # url_base
    out = ncbi_blast(sequence=df.sequence,test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["url_base"] == "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

    out = ncbi_blast(sequence=df.sequence,url_base="not_really_pg",test_run=True)
    blast_kwargs = out[0]
    assert blast_kwargs["url_base"] == "not_really_pg"

    bad_pg = [1,False,1.0,-1,None,str,"",{}]
    for b in bad_pg:
        print("passing bad url_base:",b)
        with pytest.raises(ValueError):
            out = ncbi_blast(sequence=df.sequence,url_base=b,test_run=True)

    # max_query_length, making sure sequence bits are processed correctly when
    # it's included
    out = ncbi_blast(sequence=df.sequence,test_run=True)
    out = ncbi_blast(sequence=df.sequence,max_query_length=10000,test_run=True)
    bad_int = [0,False,[],-1,None,str,"",{},int]
    for b in bad_int:
        print("passing bad max_query_length:",b)
        with pytest.raises(ValueError):
            out = ncbi_blast(sequence=df.sequence,max_query_length=b,test_run=True)

    # Get expected and actual sequence length for this df.sequence compiled
    # into >countX\nSEQUENCE\n ... format. Assumes there are less than 10
    # seqs in test dataset. Also assumes that the test sequences all have the
    # same length (160 amino acids)
    lens = [len(s) for s in df.sequence]
    expected_length = np.sum(lens) + len(lens)*len(">countX\n\n") - 1
    out = ncbi_blast(sequence=df.sequence,test_run=True)
    blast_kwargs = out[0]
    assert len(blast_kwargs["sequence"]) == expected_length

    # If our max query length is shorter than one of our sequences, it should
    # catch and die
    long_indiv_sequence = np.max([len(s) for s in df.sequence])
    with pytest.raises(ValueError):
        out = ncbi_blast(sequence=df.sequence,
                         max_query_length=long_indiv_sequence//2,
                         test_run=True)

    # Make sure splitting looks reasonable -- each sequence on own
    out = ncbi_blast(sequence=df.sequence,max_query_length=200,test_run=True)
    assert len(out) == 5
    for i,o in enumerate(out):
        assert o["sequence"].startswith(f">count{i}\n")

    # Make sure splitting looks reasonable -- 2, 2, 1
    out = ncbi_blast(sequence=df.sequence,max_query_length=340,test_run=True)
    assert len(out) == 3
    nums = [0,2,4]
    for i,o in enumerate(out):
        assert o["sequence"].startswith(f">count{nums[i]}\n")

    # Make sure splitting looks reasonable -- 2, 2, 1
    out = ncbi_blast(sequence=df.sequence,max_query_length=510,test_run=True)
    assert len(out) == 2
    nums = [0,3]
    for i,o in enumerate(out):
        assert o["sequence"].startswith(f">count{nums[i]}\n")

    # Make sure splitting looks reasonable -- none
    out = ncbi_blast(sequence=df.sequence,max_query_length=2000,test_run=True)
    assert len(out) == 1
    nums = [0]
    for i,o in enumerate(out):
        assert o["sequence"].startswith(f">count{nums[i]}\n")

    # Extra kwargs
    out = ncbi_blast(sequence=df.sequence,taxid="9606",test_run=True,
                     extra_kwarg=7)
    blast_kwargs = out[0]
    assert blast_kwargs["extra_kwarg"] == 7
