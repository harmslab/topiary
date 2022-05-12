__decription__ = \
"""
Light wrapper with muscle, compatible with muscle 3.8 and 5.1.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-05-02"

import topiary
import pandas as pd
import subprocess, sys, os, random, string

def run_muscle(input,
               output_fasta,
               muscle_binary="muscle",
               muscle_cmd_args=[]):
    """
    Run muscle for sequence alignment.

    Parameters
    ----------
    input: input to align (fasta file or topiary df)
    output_fasta: output fasta file to store alignment
    muscle_binary: location of muscle binary (default assumes 'muscle'
                   command is in the PATH).
    muscle_args: list of arguments to pass directly to muscle. Wrapper
                 specifies "-align" and "-output" (or -in/-out for old
                 version of the command line), but leaves rest as default.
                 Format should be something like ['-replicates',20,...].
                 Arguments are not checked by this function, but passed
                 directly to muscle.

    Output
    ------
    an aligned version of sequences from input in output_file. If input is
    a topiary dataframe, return a copy of the dataframe with the aligned
    sequences loaded into the alignment column.
    """

    if type(output_fasta) is not str:
        err = "\noutput_fasta '{output_fasta}' should be a string indicating\n"
        err += "fasta file where alignment should be written out.\n"
        raise ValueError(err)

    if type(input) is str:

        if not os.path.exists(input):
            err = f"\n'{input}' file does not exist.\n\n"
            raise FileNotFoundError(err)

        _run_muscle(input,output_fasta,muscle_binary)

        return None

    elif type(input) is pd.DataFrame:

        df = topiary.util.check_topiary_dataframe(input)

        # Create temporary input file
        tmp_file_root = "".join([random.choice(string.ascii_letters) for i in range(10)])
        input_fasta = "topiary-tmp_{}_blast-in.fasta".format(tmp_file_root)
        topiary.write_fasta(df,input_fasta)

        print(f"Wrote out temporary file {input_fasta} for alignment.")

        _run_muscle(input_fasta,output_fasta,muscle_binary)

        # Delete temporary input file
        try:
            os.remove(input_fasta)
        except FileNotFoundError:
            pass

        df = topiary.read_fasta_into(df,output_fasta)

        return df

    else:
        err = "\ninput should either be a fasta file or a topiary dataframe\n\n"
        raise ValueError(err)


def _run_muscle(input_fasta,
                output_fasta,
                muscle_binary="muscle"):
    """
    Run muscle for sequence alignment.

    Parameters
    ----------
        input_fasta: input fasta file to align
        output_fasta: output fasta file to store alignment
        muscle_binary: location of muscle binary (default assumes 'muscle'
                      command is in the PATH).

    Output
    ------
        an aligned version of sequences from input_file in output_file
    """

    # Make sure that muscle is in the path
    try:
        subprocess.run([muscle_binary])
    except FileNotFoundError:
        err = f"\nmuscle binary '{muscle_binary}' not found in path\n\n"
        raise ValueError(err)

    # This checks for the version of the muscle command line in use. The old
    # version used -version and will throw an error with this call; the current
    # version uses --version and this will work. This is not the most robust
    # check of versions, but is likely more robust than some kind of regex given
    # variety of compiled versions out there...
    ret_version = subprocess.run([muscle_binary,"--version"])
    if ret_version.returncode == 0:
        cmd = [muscle_binary,"-align",input_fasta,"-output",output_fasta]
    else:
        cmd = [muscle_binary,"-in",input_fasta,"-out",output_fasta]

    # This relatively complicated code is effectively subprocess.run, but
    # captures and dumps output to stdout as it is generated.
    popen = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             universal_newlines=True)
    for line in popen.stdout:
        print(line,end="")
        sys.stdout.flush()
    popen.stdout.close()
    return_code = popen.wait()

    # If running muscle failed...
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, cmd)

    print(f"\nSuccess. Output written to '{output_fasta}'.")
