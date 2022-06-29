"""
Animation indicating topiary is waiting for something without knowing
how long it will take. Launches animation on its own thread.
"""

import multiprocessing as mp
import time

class WaitingAnimation:
    """
    Animation indicating topiary is waiting for something without knowing
    how long it will take. Launches animation on its own thread.

    Parameters:
    delay : float
        time to wait between animation steps (seconds)
    num_stack : int
        how many icons to stack for one complete cycle of the animation
    icon : str
        string indicating what icon to use (can be unicode)
    """

    def __init__(self,
                 delay=1,
                 num_stack=5,
                 icon="\U0001F50D"):

        self._delay = delay
        self._num_stack = num_stack
        self._status = []
        for i in range(self._num_stack + 1):
            this_status = [" " for _ in range(self._num_stack)]
            for j in range(i):
                this_status[j] = icon
            self._status.append(" ".join(this_status))
        self._clear = num_stack*10*" "

        self._queue = mp.Queue()
        self._proc = None

    def _iterate(self,queue):

        counter = 0
        while True:

            if not queue.empty():
                break

            print(self._status[counter],end="\r",flush=True)
            time.sleep(self._delay)

            counter += 1
            if counter > self._num_stack:
                counter = 0
                print(self._clear,end="\r",flush=True)

    def start(self):
        """
        Start the animation on its own thread.
        """

        print("")
        print(self._clear,end="\r",flush=True)

        self._proc = mp.Process(target=self._iterate,args=(self._queue,))
        self._proc.start()

    def stop(self):
        """
        Stop the animation.
        """

        self._queue.put(True)
        self._proc.join()

        print(self._clear,end="\r",flush=True)
        print("",flush=True)
