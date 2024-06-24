
import pytest

from topiary._private import threads
from topiary._private.threads import _thread
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

def _function_to_thread(value,some_option=False,lock=None):
    """
    Function for testing thread_manager.
    """

    if lock is not None:

        # If lock is acquired, the output value will change. (Not directly
        # testing locking -- that will depend too much on specific function).
        lock.acquire()
        try:
            value = 2*value
        finally:
            lock.release()

    if some_option:
        value = 3*value

    # Sleep for a random time to knock multithreaded-order out of sync
    time.sleep(random.choice([0,0.01]))

    return value

def test_thread_manager():

    kwargs_list = [{"value":2,"some_option":False},
                   {"value":2,"some_option":True},
                   {"value":3,"some_option":False}]

    # Run on single thread. (This will not use pool)
    output = threads.thread_manager(kwargs_list,
                                    _function_to_thread,
                                    num_threads=1,
                                    progress_bar=False)
    assert np.array_equal(output,[2,6,3])

    # Run on single thread. (This will not use pool; make sure progress bar at
    # least accepted).
    output = threads.thread_manager(kwargs_list,
                                    _function_to_thread,
                                    num_threads=1,
                                    progress_bar=True)
    assert np.array_equal(output,[2,6,3])

    # Run twice on two threads to make sure order is preserved
    for i in range(2):
        output = threads.thread_manager(kwargs_list,
                                        _function_to_thread,
                                        num_threads=2,
                                        progress_bar=True,
                                        pass_lock=False)
        assert np.array_equal(output,[2,6,3])

    # Run twice on two threads make sure order is preserved
    for i in range(2):
        output = threads.thread_manager(kwargs_list,
                                        _function_to_thread,
                                        num_threads=2,
                                        progress_bar=False,
                                        pass_lock=True)
        assert np.array_equal(output,[4,12,6])


def test__thread():
    # tested implicitly in test_thread_manager. touch but do not call so test
    # crawler does not flag as missing
    _thread
    return None
