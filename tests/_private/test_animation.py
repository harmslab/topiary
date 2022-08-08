
import pytest

from topiary._private import animation

import time

def test_WaitingAnimation():

    # Not really testing animation output, but am testing threading. This test
    # really only makes sense if we're using pytest-timeout with the
    # --timeout flag. If threading doesn't start properly, time.sleep will not
    # occur; of w.stop() does not exit, this function will hang forever. 

    w = animation.WaitingAnimation()

    w.start()
    time.sleep(2)
    w.stop()
def test_WaitingAnimation__iterate():

    pass

def test_WaitingAnimation_start():

    pass

def test_WaitingAnimation_stop():

    pass

