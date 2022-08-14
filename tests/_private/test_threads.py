
import pytest

from topiary._private import threads
import time, random

import numpy as np

def test_MockLock():

    pass

def test_get_num_threads():

    num_threads = threads.get_num_threads(-1,manual_num_cores=10)
    assert num_threads == 10

    num_threads = threads.get_num_threads(1,manual_num_cores=10)
    assert num_threads == 1

    num_threads = threads.get_num_threads(5,manual_num_cores=10)
    assert num_threads == 5

    num_threads = threads.get_num_threads(5,manual_num_cores=2)
    assert num_threads == 2

    bad_num_threads = [0,-2,"test",None,int,1.1]
    for b in bad_num_threads:
        print(f"testing {b} in get_num_threads")
        with pytest.raises(ValueError):
            threads.get_num_threads(b,manual_num_cores=10)

def _function_to_thread(rocking,usa=False,lock=None):
    """
    Function for testing thread_manager.
    """

    if lock is not None:

        # If lock is acquired, the output value will change. (Not directly
        # testing locking -- that will depend too much on specific function).
        lock.acquire()
        try:
            rocking = 2*rocking
        finally:
            lock.release()

    if usa:
        rocking = 3*rocking

    # Sleep for a random time to knock multithreaded-order out of sync
    time.sleep(random.choice([0,0.2,0.4,0.6,0.8]))

    return rocking

def test_thread_manager():

    kwargs_list = [{"rocking":2,"usa":False},
                   {"rocking":2,"usa":True},
                   {"rocking":3,"usa":False},
                   {"rocking":3,"usa":True}]

    # Run on single thread. (This will not use pool)
    output = threads.thread_manager(kwargs_list,
                                    _function_to_thread,
                                    num_threads=1,
                                    progress_bar=False)
    assert np.array_equal(output,[2,6,3,9])

    # Run on single thread. (This will not use pool; make sure progress bar at
    # least accepted).
    output = threads.thread_manager(kwargs_list,
                                    _function_to_thread,
                                    num_threads=1,
                                    progress_bar=True)
    assert np.array_equal(output,[2,6,3,9])

    # Run five times on two threads to make sure order is preserved
    for i in range(5):
        output = threads.thread_manager(kwargs_list,
                                        _function_to_thread,
                                        num_threads=2,
                                        progress_bar=True,
                                        pass_lock=False)
        assert np.array_equal(output,[2,6,3,9])

    # Run five times to make sure order is preserved
    for i in range(5):
        output = threads.thread_manager(kwargs_list,
                                        _function_to_thread,
                                        num_threads=2,
                                        progress_bar=False,
                                        pass_lock=True)
        assert np.array_equal(output,[4,12,6,18])

def test__thread():
    # tested implicitly in test_thread_manager.
    return True
