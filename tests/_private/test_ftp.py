
import pytest

from topiary._private.ftp import _ftp_thread
from topiary._private.ftp import ftp_download
from topiary._private.ftp import calc_md5

import os
import sys

def test__ftp_thread(tmpdir):

    cwd = os.getcwd()
    os.chdir(tmpdir)

    _ftp_thread(file_name="GCF_000001405.40_GRCh38.p14_assembly_report.txt",
                path="/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14",
                url="ftp.ncbi.nlm.nih.gov",
                kwargs={})
    assert os.path.isfile("GCF_000001405.40_GRCh38.p14_assembly_report.txt")

    os.chdir(cwd)

def test_ftp_download(tmpdir):

    def _get_creation_time(file_name):

        if sys.platform == "darwin":
            creation_time = os.stat(file_name).st_birthtime
        else:
            creation_time = None

        return creation_time

    cwd = os.getcwd()
    os.chdir(tmpdir)

    # Download

    file_name = "GCF_000001405.40_GRCh38.p14_assembly_report.txt"

    ftp_download(file_name=file_name,
                 path="/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14",
                 url="ftp.ncbi.nlm.nih.gov")
    assert os.path.isfile("GCF_000001405.40_GRCh38.p14_assembly_report.txt")

    creation_time = _get_creation_time(file_name)

    # Should just keep existing

    ftp_download(file_name=file_name,
                 path="/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14",
                 url="ftp.ncbi.nlm.nih.gov",
                 resume=True)
    assert os.path.isfile("GCF_000001405.40_GRCh38.p14_assembly_report.txt")

    new_creation_time = _get_creation_time(file_name)
    if creation_time is not None:
        assert new_creation_time == creation_time

    # Should overwrite existing

    ftp_download(file_name=file_name,
                 path="/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14",
                 url="ftp.ncbi.nlm.nih.gov",
                 resume=False)
    assert os.path.isfile("GCF_000001405.40_GRCh38.p14_assembly_report.txt")

    new_creation_time = _get_creation_time(file_name)
    if creation_time is not None:
        assert new_creation_time != creation_time

    os.chdir(cwd)

def test_calc_md5(ftp_test_files,tmpdir):

    assert calc_md5(ftp_test_files["GCF_000001405.40_GRCh38.p14_assembly_report.txt.gz"]) == "edbd21ee24986ce383fb54d2a7f93708"
