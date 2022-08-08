"""
Functions for interacting with ftp servers.
"""

import numpy as np

from tqdm.auto import tqdm

import ftplib
import os
import hashlib
import time
import multiprocessing as mp


def _ftp_thread(file_name,path,url,kwargs):
    """
    Download a file off an ftp server on its own thread.

    Parameters
    ----------
    file_name : str
        name of file on server
    path : str
        path to file on server
    url : str
        url pointing to server
    kwargs : dict
        kwargs to send to ftp.retrbinary
    """

    ftp = ftplib.FTP(url)
    ftp.login("anonymous", "")
    ftp.cwd(path)

    ftp.retrbinary(cmd="RETR " + file_name,
                   callback=open(file_name, 'ab').write,
                   **kwargs)

    ftp.quit()

def ftp_download(file_name,path,url,resume=True,silent=False):
    """
    Download a file by anonymous ftp, resuming halted download if possible.

    Parameters
    ----------
    file_name : str
        name of file on server
    path : str
        path to file on server
    url : str
        url pointing to server
    resume : bool, default=True
        resume halted download if possible
    silent : bool, default=False
        whether or not to print status and progress bar
    """

    if not silent:
        print(f"Downloading {file_name}",flush=True)

    local_file_size = 0
    kwargs = {}
    if os.path.exists(file_name):
        if resume:
            print("Local file found. Checking to resume partial download.",flush=True)
            local_file_size = os.path.getsize(file_name)
            kwargs["rest"] = str(local_file_size)
        else:
            os.remove(file_name)

    # Get the size of the file off the server
    ftp = ftplib.FTP(url)
    ftp.login("anonymous", "")
    ftp.cwd(path)
    remote_file_size = ftp.size(file_name)
    ftp.quit()

    # Launch download on a thread
    download_proc = mp.Process(target=_ftp_thread,args=(file_name,path,url,kwargs))
    download_proc.start()

    if not silent:

        # Continuously monitor the size of the file and update progress bar
        size_in_mb = np.round(remote_file_size/1e6,1)
        with tqdm(total=size_in_mb) as pbar:

            download_done = False
            while not download_done:
                if os.path.exists(file_name):
                    current_file_size = np.round(os.path.getsize(file_name)/1e6,1)
                    if current_file_size >= size_in_mb:
                        current_file_size = size_in_mb
                        download_done = True

                    pbar.n = current_file_size
                    pbar.refresh()

                time.sleep(0.1)

    # Join downlod thread
    download_proc.join()


def calc_md5(file_name,chunk_size=4096):
    """
    Get the md5 hash string for a file.

    Parameters
    ----------
    file_name : str
        file on which to calculate md5
    chunk_size : int, default=4096
        read file in chunks of this size

    Returns
    -------
    md5 : str
        md5 hexdigest for file
    """

    hash_md5 = hashlib.md5()
    with open(file_name, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
