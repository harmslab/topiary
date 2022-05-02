
import pytest
from topiary.ncbi import ncbi_blast

import numpy as np

def test_ncbi_blast(test_dataframes):

    df = test_dataframes["good-df"]
    blast_kwargs = ncbi_blast(sequence=df.sequence,test_run=True)

    # Make sure it got the sequences right
    sequences = blast_kwargs["sequence"]
    sequences = sequences.split("\n")
    sequences = [sequences[i] for i in range(1,len(sequences),2)]

    assert np.array_equal(sequences,df.sequence)

    # Make sure other arguments are processing properly
    blast_kwargs = ncbi_blast(sequence=df.sequence,test_run=True)
    assert blast_kwargs["database"] == "nr"
    assert blast_kwargs["hitlist_size"] == '50'
    assert blast_kwargs["program"] == "blastp"
    assert blast_kwargs["expect"] == '0.01'
    assert blast_kwargs["gapcosts"] == '11 1'
    assert blast_kwargs["url_base"] == "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    with pytest.raises(KeyError):
        blast_kwargs["entrez_query"]

    # test sequence bits
    blast_kwargs = ncbi_blast(sequence="test",test_run=True)
    sequences = blast_kwargs["sequence"]
    assert sequences.split("\n")[1] == "test"

    blast_kwargs = ncbi_blast(sequence=["test"],test_run=True)
    sequences = blast_kwargs["sequence"]
    assert sequences.split("\n")[1] == "test"

    blast_kwargs = ncbi_blast(sequence=["test","this"],test_run=True)
    sequences = blast_kwargs["sequence"]
    assert sequences.split("\n")[1] == "test"
    assert sequences.split("\n")[3] == "this"

    # Send in bad sequences
    bad_sequences = [1,False,1.0,-1,None,str,""]
    for b in bad_sequences:
        print("passing bad sequence:",b)
        with pytest.raises(ValueError):
            blast_kwargs = ncbi_blast(sequence=b,test_run=True)
        with pytest.raises(ValueError):
            blast_kwargs = ncbi_blast(sequence=[b],test_run=True)

    # taxid
    blast_kwargs = ncbi_blast(sequence=df.sequence,taxid=9606,test_run=True)
    assert blast_kwargs["entrez_query"] == "txid9606[ORGN]"

    blast_kwargs = ncbi_blast(sequence=df.sequence,taxid="9606",test_run=True)
    assert blast_kwargs["entrez_query"] == "txid9606[ORGN]"

    blast_kwargs = ncbi_blast(sequence=df.sequence,taxid=(9606,1234),test_run=True)
    assert blast_kwargs["entrez_query"] == "txid9606[ORGN] or txid1234[ORGN]"

    blast_kwargs = ncbi_blast(sequence=df.sequence,taxid=["9606","1234"],test_run=True)
    assert blast_kwargs["entrez_query"] == "txid9606[ORGN] or txid1234[ORGN]"

    bad_taxid = [1.15,str,int]
    for b in bad_taxid:
        with pytest.raises(ValueError):
            print("passing bad taxid:",b)
            blast_kwargs = ncbi_blast(sequence=df.sequence,taxid=b,test_run=True)

    # Final check of kwargs
    blast_kwargs = ncbi_blast(sequence=df.sequence,taxid="9606",test_run=True,
                              extra_kwarg=7)
    assert blast_kwargs["extra_kwarg"] == 7
