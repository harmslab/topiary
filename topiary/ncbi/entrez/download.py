"""
Functions to download files off of the NCBI via FTP.
"""

from topiary._private import check
from topiary._private.ftp import ftp_download
from topiary._private.ftp import calc_md5

import os

def _read_md5_file(md5_file):
    """
    Read an NCBI md5 checksum file.

    Parameters
    ----------
    md5_file : str
        file with md5sum entries

    Returns
    -------
    md5_dict : dict
        dictionary keying file names to md5 hex strings
    """

    md5_dict = {}
    with open(md5_file) as f:
        for line in f:
            if line.strip() == "":
                continue

            col = line.split()
            file = col[1][2:].strip()

            md5_dict[file] = col[0].strip()

    return md5_dict


def ncbi_ftp_download(full_url,
                      file_base="_protein.faa.gz",
                      md5_file="md5checksums.txt",
                      num_attempts=5):
    """
    Download a proteome, genome, etc. from the ncbi. Makes multiple attempts,
    resumes interrupted downloads, and checks md5sum to validate integrity.

    Parameters
    ----------
    full_url : str
        full url to directory containing data. This will look like
        ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14
        (leading ftp:// allowed).
    file_base : str, default="_protein.faa.gz"
        class of file to download. "_protein.faa.gz" corresponds to a proteome,
        "_genomic.fna.gz" to a genome, etc.
    md5_file : str, default="md5checksums.txt"
        name of md5 checksum file on the server.
    num_attempts : int, default=5
        number of times to try to download before giving up
    """

    full_url = str(full_url)
    if full_url.startswith("ftp://"):
        full_url = full_url[6:]

    file_base = str(file_base)
    md5_file = str(md5_file)
    num_attempts = check.check_int(num_attempts,
                                   "num_attempts",
                                   minimum_allowed=1)

    split = full_url.split("/")
    url = split[0]
    path = "/" + "/".join([s for s in split[1:] if s != ""])
    file_name = f"{split[-1]}{file_base}"

    # Get md5 sum file
    md5_attempts = 0
    while True:
        try:
            ftp_download(md5_file,path,url,resume=False,silent=True)
            md5_dict = _read_md5_file(md5_file)
            break
        except Exception as e:
            md5_attempts += 1
            if md5_attempts > 5:
                err = "Could not download md5sum file.\n"
                raise Exception(err) from e

    try:
        md5_dict[file_name]
    except KeyError:
        err = f"The file '{file_name}' is not present on the NCBI.\n"
        err += f"Full path: {full_url}/{path}/{file_name}"
        raise FileNotFoundError(err)

    counter = 0
    success = False
    while counter < num_attempts:

        print(f"Attempt {counter + 1} of {num_attempts}.")
        ftp_download(file_name,path,url,resume=True,silent=False)

        # Make sure it downloaded
        if not os.path.isfile(file_name):
            counter += 1
            continue

        # Check md5sum
        file_md5 = calc_md5(file_name)
        if md5_dict[file_name] != file_md5:
            print("md5 hash does not match.")
            print(f"Expected {md5_dict[file_name]}")
            print(f"Got {file_md5}")
            print("Trying again...")

            os.remove(file_name)

            counter += 1
            continue

        success = True
        break

    if not success:
        err = f"\n\nCould not download {url}/{path}/{file_name}.\n"
        err += "This is likely due to server congestion and/or download limits.\n"
        err += "Please try again later.\n\n"
        raise RuntimeError(err)
