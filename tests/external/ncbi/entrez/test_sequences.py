import pytest
import topiary

from topiary.external.ncbi.entrez.sequences import _get_sequences_thread_function
from topiary.external.ncbi.entrez.sequences import get_sequences

def test__get_sequences_thread_function():

    import multiprocessing as mp
    lock = mp.Lock()
    results = _get_sequences_thread_function(ids="EAW87024.1,PJI62056",
                                             num_tries_allowed=5,
                                             lock=lock)

    assert isinstance(results,list)
    assert len(results) == 2
    assert results[0][0].split(".")[0] == "EAW87024"
    assert results[0][1][:10] == "MLPFLFFSTL"
    assert results[1][0].split(".")[0] == "PJI62056"
    assert results[1][1][:10] == "MKPVTLYDVA"

def test_get_sequences():

    to_download = ["EAW87024.1","PJI62056"]

    # Download these sequences
    results = get_sequences(to_download,block_size=50)

    assert isinstance(results,list)
    assert len(results) == 2
    assert results[0][0].split(".")[0] == "EAW87024"
    assert results[0][1][:10] == "MLPFLFFSTL"
    assert results[1][0].split(".")[0] == "PJI62056"
    assert results[1][1][:10] == "MKPVTLYDVA"

    # Force get_sequences to do queries and merge.
    results = get_sequences(to_download,block_size=1)

    assert isinstance(results,list)
    assert len(results) == 2
    assert results[0][0].split(".")[0] == "EAW87024"
    assert results[0][1][:10] == "MLPFLFFSTL"
    assert results[1][0].split(".")[0] == "PJI62056"
    assert results[1][1][:10] == "MKPVTLYDVA"

    results = get_sequences([],block_size=1)
    assert isinstance(results,list)
    assert len(results) == 0

    bad_args = [1,1.1,None,list]
    for b in bad_args:
        print("Trying",b)
        with pytest.raises(ValueError):
            get_sequences(to_download=b)


    bad_args = [-1,0,int,[1,2,3],1.1,None]
    for b in bad_args:
        print("Trying",b)

        with pytest.raises(ValueError):
            get_sequences(to_download,block_size=b)

        with pytest.raises(ValueError):
            get_sequences(to_download,num_tries_allowed=b)

        # This should work for num_threads, so skip
        if b == -1:
            continue

        with pytest.raises(ValueError):
            get_sequences(to_download,num_threads=b)
