import pytest

from topiary.ncbi.entrez.download import _read_md5_file
from topiary.ncbi.entrez.download import ncbi_ftp_download
from topiary._private.ftp import calc_md5

import os

def test__read_md5_file(ftp_test_files):

    # Make sure we can read an md5checksums file

    md5test_file = ftp_test_files["md5checksums_test.txt"]

    expected_length = 0
    with open(md5test_file) as f:
        for line in f:
            expected_length += 1

    out_dict = _read_md5_file(md5test_file)
    assert isinstance(out_dict,dict)

    assert out_dict["README_patch_release.txt"] == "3b7f12ebd3d129698e86fb8701bb9688"
    assert out_dict["GCF_000001405.40_GRCh38.p14_protein.faa.gz"] == "6dbdf8c7f9c7f39da15abb996831c733"
    assert len(out_dict) == 1419

def test_ncbi_ftp_download(ftp_test_files,tmpdir):

    cwd = os.getcwd()
    os.chdir(tmpdir)

    test_url = ["ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14",
                "ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14"]
    file_base = "_assembly_report.txt"

    for u in test_url:

        ncbi_ftp_download(u,
                          file_base=file_base,
                          md5_file="md5checksums.txt",
                          num_attempts=5)
        assert os.path.isfile("md5checksums.txt")
        downloaded = f"GCF_000001405.40_GRCh38.p14{file_base}"
        assert os.path.isfile(downloaded)
        md5_dict = _read_md5_file("md5checksums.txt")
        assert calc_md5(downloaded) == md5_dict[downloaded]


    with pytest.raises(ValueError):
        ncbi_ftp_download(test_url[0],
                          file_base=file_base,
                          md5_file="md5checksums.txt",
                          num_attempts=-1)


    os.chdir(cwd)
