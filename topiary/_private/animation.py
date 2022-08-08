"""
Animation indicating topiary is waiting for something without knowing
how long it will take. Launches animation on its own thread.
"""

import topiary

import multiprocessing as mp
import time, sys

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
            this_status.append("\r")
            self._status.append(" ".join(this_status))
        self._clear = num_stack*10*" " + "\r"

        self._stop_queue = mp.Queue()
        self._proc = None

    def _iterate(self,stop_queue):

        counter = 0
        while True:

            sys.stdout.write(self._status[counter])
            sys.stdout.flush()
            time.sleep(self._delay)

            if not stop_queue.empty():
                break

            counter += 1
            if counter > self._num_stack:
                counter = 0
                sys.stdout.write(self._clear)

    def start(self):
        """
        Start the animation on its own thread.
        """

        sys.stdout.write("\n")
        sys.stdout.write(self._clear)
        sys.stdout.flush()

        self._proc = mp.Process(target=self._iterate,args=(self._stop_queue,))
        self._proc.start()

    def stop(self):
        """
        Stop the animation.
        """

        self._stop_queue.put(True)
        self._proc.join()

        sys.stdout.write(self._clear)
        sys.stdout.write("\n")
        sys.stdout.flush()
