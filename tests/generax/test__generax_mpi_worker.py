
import pytest

# DO NOT IMPORT THESE MODULES. If these are imported, mpirun just doesn't do...
# anything.
#import topiary.generax._generax_mpi_worker
#from topiary.generax._generax_mpi_worker import main

import shutil
import subprocess
import os
import glob
import sys

SCRIPT = os.path.join(os.path.dirname(__file__),
                      "..","..","topiary","generax","_generax_mpi_worker.py")

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_main(generax_data,tmpdir):

    def _check_dir(dir_name,num_workers,run_id,already_complete=[]):


        workers_seen = []
        for d in os.listdir(dir_name):

            this_dir = os.path.join(dir_name,d)
            if not os.path.isdir(this_dir):
                continue

            if d in already_complete:

                # Other completed file - not one we made!
                completed_files = glob.glob(os.path.join(this_dir,"completed*"))
                assert len(completed_files) == 1
                assert not completed_files[0].startswith(f"completed_run-{run_id}")

                # Should not have run because already complted
                assert not os.path.isfile(os.path.join(this_dir,
                                                       "generax-ran.txt"))

                continue

            # Only one complted file
            completed_files = glob.glob(os.path.join(this_dir,"completed*"))
            assert len(completed_files) == 1

            # File has expected structure
            parts = os.path.split(completed_files[0])[-1].split("_")
            assert parts[0] == "completed"
            assert parts[1] == f"run-{run_id}"
            assert parts[2].startswith("by-")

            # Make sure worker is correct
            worker = int(parts[2].split("-")[1])
            assert worker <= (num_workers - 1)
            workers_seen.append(worker)

            # Make sure thing rang
            assert os.path.isfile(os.path.join(this_dir,"generax-ran.txt"))

        # Make sure all workers ran
        assert len(set(workers_seen)) == num_workers

    # Get path to mpirun so we know what it is if bad things happen
    mpirun = shutil.which("mpirun")

    test_dir = generax_data["test-dir"]
    gen_tmp = os.path.join(tmpdir,"generax")
    os.mkdir(gen_tmp)

    for i in [1,2]:

        this_test_dir = os.path.join(gen_tmp,f"test_worker_{i}")
        shutil.copytree(test_dir,this_test_dir)

        run_id = f"testid{i}"
        cmd = [mpirun,"-np",f"{i}",sys.executable,SCRIPT,this_test_dir,run_id]

        print("Running:"," ".join(cmd))
        ret = subprocess.run(cmd,capture_output=True)

        _check_dir(this_test_dir,i,run_id)

    for i in [1,2]:

        this_test_dir = os.path.join(gen_tmp,f"test_worker_{i}_inject_claim")
        shutil.copytree(test_dir,this_test_dir)

        # Inject a competing run file but *not* complete into directory 0.
        # Everything should run, leaving this file alone.
        other_claim = os.path.join(this_test_dir,"0","claimed_run-other_20")
        f = open(other_claim,"w")
        f.write("")
        f.close()

        run_id = f"testid{i}"
        cmd = [mpirun,"-np",f"{i}",sys.executable,SCRIPT,this_test_dir,run_id]

        print("Running:"," ".join(cmd))
        ret = subprocess.run(cmd,capture_output=True)

        _check_dir(this_test_dir,i,run_id)

        assert os.path.isfile(other_claim)

    for i in [1,2]:

        this_test_dir = os.path.join(gen_tmp,f"test_worker_{i}_inject_completed")
        shutil.copytree(test_dir,this_test_dir)

        # Inject a competing complete. Should run others but leave this alone.
        other_completed = os.path.join(this_test_dir,"0","completed")
        f = open(other_completed,"w")
        f.write("")
        f.close()

        run_id = f"testid{i}"
        cmd = [mpirun,"-np",f"{i}",sys.executable,SCRIPT,this_test_dir,run_id]

        print("Running:"," ".join(cmd))
        ret = subprocess.run(cmd,capture_output=True)

        _check_dir(this_test_dir,i,run_id,already_complete=["0"])

        # Should not have touched this other complete file
        assert os.path.isfile(other_completed)
