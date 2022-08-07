
import pytest

import topiary
from topiary.external.ncbi.blast.make import make_blast_db

import os
import glob

def test_make_blast_db(make_blast_db_files,tmpdir):

    # Skip test on windows machine
    if os.name == "nt":
        print("Skipping creation of database test on windows.")
        return

    cwd = os.getcwd()
    os.chdir(tmpdir)
    expected_extensions = ["pin","phr","psq"]

    faa = [make_blast_db_files["test1.faa"],make_blast_db_files["test2.faa"]]
    out = "output"
    make_blast_db(faa,db_name=out)
    print("Check:")
    print(os.listdir("."))
    for e in expected_extensions:
        out_file = f"{out}.{e}"
        assert os.path.isfile(out_file)

    for out_file in glob.glob("output.p*"):
        os.remove(out_file)

    faa_gz = [make_blast_db_files["test1.faa.gz"],make_blast_db_files["test2.faa.gz"]]
    make_blast_db(faa_gz,db_name=out)
    print("Check:")
    print(os.listdir("."))
    for e in expected_extensions:
        out_file = f"{out}.{e}"
        assert os.path.isfile(out_file)

    for out_file in glob.glob("output.p*"):
        os.remove(out_file)

    mixed = [make_blast_db_files["test1.faa.gz"],make_blast_db_files["test2.faa"]]
    make_blast_db(mixed,db_name=out)
    print("Check:")
    print(os.listdir("."))
    for e in expected_extensions:
        out_file = f"{out}.{e}"
        assert os.path.isfile(out_file)

    # NOTE: did *not* remove files from previous test.
    with pytest.raises(FileExistsError):
        make_blast_db(faa,db_name=out)

    make_blast_db(faa,db_name=out,overwrite=True)

    # Make sure it's really a blast database...
    topiary.ncbi.blast.local_blast(sequence="SOMESEQVENCE",db=out)

    os.chdir(cwd)
