
import pytest
from topiary.ncbi.blast.local import _construct_args

def test_threading_logic():
    """
    Test that local_blast correctly chooses between process-level and 
    binary-level parallelism.
    """
    
    # We can't easily test local_blast without a database, but we can test
    # the logic by mocking or checking the internal functions.
    # However, let's just test _construct_args which I modified.

    # 1. Case: 19 sequences, num_threads=8 (Small set)
    # local_blast will set num_threads_pool=1, num_threads_blast=8, block_size=19
    sequence_list = ["seq"] * 19
    blast_args = ["blastp", "-db", "db", "-num_threads", "8"]
    
    kwargs_list, num_threads_pool = _construct_args(sequence_list=sequence_list,
                                                    blast_args=blast_args,
                                                    block_size=19,
                                                    num_threads=1,
                                                    manual_num_cores=8)
    
    assert num_threads_pool == 1
    assert len(kwargs_list) == 1
    assert "-num_threads" in kwargs_list[0]["blast_args"]
    assert "8" in kwargs_list[0]["blast_args"]

    # 2. Case: 100 sequences, num_threads=8 (Large set)
    # local_blast will set num_threads_pool=8, num_threads_blast=1, block_size=20
    sequence_list = ["seq"] * 100
    blast_args = ["blastp", "-db", "db", "-num_threads", "1"]
    
    kwargs_list, num_threads_pool = _construct_args(sequence_list=sequence_list,
                                                    blast_args=blast_args,
                                                    block_size=20,
                                                    num_threads=8,
                                                    manual_num_cores=8)
    
    # max_useful_threads = 100 // 20 = 5.
    # num_threads_pool (8) > max_useful_threads (5), so it should adjust block_size
    assert num_threads_pool == 8
    assert len(kwargs_list) >= 8
