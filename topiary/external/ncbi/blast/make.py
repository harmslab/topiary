"""
Construct a blast database.
"""

import topiary
from topiary._private import check

import os, subprocess, gzip, random, string

def make_blast_db(input_files,db_name,overwrite=False,makeblastdb_binary="makeblastdb"):
    """
    Make a protein blast database from a set of fasta input files.

    Parameters
    ----------
    input_files : list
        list of input files. Must be fasta files (with extension .faa) or
        zipped fasta files (with extension .faa.gz).
    db_name : str
        name of final database
    overwrite : bool, default=False
        whether or not to overwrite an existing database
    makeblastdb_binary : str, default="makeblastdb"
        makeblastdb binary to use

    Returns
    -------
    None
    """

    print("\nCreating blast database.\n",flush=True)

    # If the user passes in a string for input_files, make it into a list
    if type(input_files) is str:
        input_files = [input_files]

    if not hasattr(input_files,"__iter__") or type(input_files) is type:
        err = "\ninput_files should be a list of input files to put into the\n"
        err += "blast database.\n\n"
        raise ValueError(err)

    # Check for sanity of db_name
    db_name = str(db_name)

    # Process overwrite argument
    overwrite = check.check_bool(overwrite)

    # If we are overwriting an existing blast database.
    if os.path.exists(f"{db_name}.psq"):
        if overwrite:
            expected_extensions = ["pdb","pin","pot","ptf","phr","psq","pto"]
            for e in expected_extensions:
                try:
                    os.remove(f"{db_name}.{e}")
                except FileNotFoundError:
                    pass
        else:
            err = f"\nA blast database with name {db_name} already exists.\n\n"
            raise FileExistsError(err)

    # Make sure that makeblastdb is in the path
    try:
        subprocess.run([makeblastdb_binary],capture_output=True)
    except FileNotFoundError:
        err = f"\makeblastdb_binary binary '{makeblastdb_binary}' not found in path\n\n"
        raise ValueError(err)

    rand_string = "".join([random.choice(string.ascii_letters) for _ in range(10)])
    tmp_faa_name = f"topiary_makeblastdb_{rand_string}.faa"

    # Combine files into a single temporary file
    out = open(tmp_faa_name,'a')
    for i, file in enumerate(input_files):
        print(f"Reading {file}",flush=True)

        if file[-7:] == ".faa.gz":
            with gzip.open(str(file)) as f:
                for line in f:
                    out.write(line.decode())
        elif file[-4:] == ".faa":
            with open(str(file)) as f:
                for line in f:
                    out.write(line)
        else:
            out.close()
            os.remove(tmp_faa_name)
            err = "\nInput files must have the extension .faa or .faa.gz.\n\n"
            raise ValueError(err)

    out.close()



    # Simple makeblastdb call
    cmd = [makeblastdb_binary,"-dbtype","prot","-in",tmp_faa_name,"-out",db_name]

    # This relatively complicated code is effectively subprocess.run, but
    # captures and dumps output to stdout as it is generated. This is useful
    # for a jupyter notebook, as the output appears in the notebook rather than
    # on the command line.
    popen = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             universal_newlines=True)
    for line in popen.stdout:
        print(line,end="",flush=True)
    popen.stdout.close()
    return_code = popen.wait()

    os.remove(tmp_faa_name)

    print("Done.",flush=True)
