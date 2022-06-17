__decription__ = \
"""
Light wrapper for muscle, compatible with muscle 3.8 and 5.1.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-05-02"

import topiary
from topiary import _arg_processors
import pandas as pd
import subprocess, sys, os, random, string, re

def run_muscle(input,
               output_fasta=None,
               super5=False,
               muscle_cmd_args=[],
               muscle_binary="muscle"):
    """
    Run muscle to align sequences.

    Parameters
    ----------
        input: input to align (fasta file or topiary df)
        output_fasta: output fasta file to store alignment. Optional if the input
                      is a dataframe; reqiured if the input is a fasta file.
        super5: bool. User the 'super5' mode of muscle 5
        muscle_cmd_args: list of arguments to pass directly to muscle. Wrapper
                         specifies "-align" and "-output" (or -in/-out for old
                         version of the command line), but leaves rest as default.
                         Format should be something like ['-replicates',20,...].
                         Arguments are not checked by this function, but passed
                         directly to muscle.
        muscle_binary: location of muscle binary (default assumes 'muscle'
                       command is in the PATH).

    Return
    ------
        If input is a topiary dataframe, return a copy of the dataframe with the
        aligned sequences loaded into the alignment column. Otherwise, write to the
        output_fasta file and return None from the function.
    """

    # If output_fasta is defined, make sure it's a string
    if output_fasta and type(output_fasta) is not str:
        err = "\noutput_fasta '{output_fasta}' should be a string indicating\n"
        err += "fasta file where alignment should be written out.\n"
        raise ValueError(err)

    # If the input is a string.
    if type(input) is str:

        # Make sure the input file exists.
        if not os.path.exists(input):
            err = f"\n'{input}' file does not exist.\n\n"
            raise FileNotFoundError(err)

        # Make sure the output file is not None.
        if output_fasta is None:
            err = "\nIf aligning a fasta file, an output_fasta file must be\n"
            err += "specified.\n"
            raise ValueError(err)

        # Do the alignment.
        _run_muscle(input,output_fasta,super5,muscle_cmd_args,muscle_binary)

        print(f"\nSuccess. Alignment written to '{output_fasta}'.",flush=True)

        return None

    # If the input is a dataframe
    elif type(input) is pd.DataFrame:

        # Check the dataframe
        df = _arg_processors.process_topiary_dataframe(input)

        # Create temporary input file
        tmp_file_root = "".join([random.choice(string.ascii_letters) for i in range(10)])
        input_fasta = "topiary-tmp_{}_align-in.fasta".format(tmp_file_root)
        topiary.write_fasta(df,input_fasta)

        # Create temporary output file name if required.
        temporary_output = False
        if output_fasta is None:
            output_fasta = "topiary-tmp_{}_align-out.fasta".format(tmp_file_root)
            temporary_output = True

        # Do the alignment
        _run_muscle(input_fasta,output_fasta,super5,muscle_cmd_args,muscle_binary)

        # Read alignment back into the dataframe
        df = topiary.read_fasta_into(df,output_fasta)

        # Delete temporary input file
        try:
            os.remove(input_fasta)
        except FileNotFoundError:
            pass

        # Delete temporary output file
        print("\nSuccess. Alignment written to the `alignment` column in the dataframe.",flush=True)
        if temporary_output:
            try:
                os.remove(output_fasta)
            except FileNotFoundError:
                pass
        else:
            print(f"\nSuccess. Alignment written to '{output_fasta}'.",flush=True)

        return df

    else:
        err = "\ninput should either be a fasta file or a topiary dataframe\n\n"
        raise ValueError(err)


def _run_muscle(input_fasta,
                output_fasta,
                super5=False,
                muscle_cmd_args=[],
                muscle_binary="muscle"):
    """
    Run muscle for sequence alignment.

    Parameters
    ----------
        input_fasta: input fasta file to align
        output_fasta: output fasta file to store alignment
        super5: bool. User the 'super5' mode of muscle 5
        muscle_cmd_args: list of arguments to pass directly to muscle. Wrapper
                         specifies "-align" and "-output" (or -in/-out for old
                         version of the command line), but leaves rest as default.
                         Format should be something like ['-replicates',20,...].
                         Arguments are not checked by this function, but passed
                         directly to muscle.
        muscle_binary: location of muscle binary (default assumes 'muscle'
                       command is in the PATH).

    Return
    ------
        None
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
        if super5:
            cmd[1] = "-super5"
    else:
        cmd = [muscle_binary,"-in",input_fasta,"-out",output_fasta]


    cmd.extend(muscle_cmd_args)

    # This relatively complicated code is effectively subprocess.run, but
    # captures and dumps output to stdout as it is generated. This is useful
    # for a jupyter notebook, as the output appears in the notebook rather than
    # on the command line.
    popen = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             universal_newlines=True)
    for line in popen.stdout:
        if re.search("100.0%",line):
            endl = "\n"
        else:
            print(90*" ",end="\r")
            endl = "\r"
        print(line.strip(),end=endl,flush=True)
    popen.stdout.close()
    return_code = popen.wait()

    # If running muscle failed...
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, cmd)
