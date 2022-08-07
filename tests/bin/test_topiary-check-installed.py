

import sys
import os
import subprocess

def test_main():

    # simple test. Make sure it runs. We test validation stack directly
    # elsewhere.

    # Get location of binary
    location = os.path.dirname(os.path.realpath(__file__))
    test_bin = os.path.join(location,"..","..","bin","topiary-check-installed")

    if os.name == "nt":
        base_cmd = [sys.executable,test_bin]
    else:
        base_cmd = [test_bin]

    ret = subprocess.run(base_cmd)
    assert ret.returncode == 0
