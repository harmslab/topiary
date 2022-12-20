"""
Interface to muscle, compatible with muscle 3.8 and 5.1.
"""

import topiary
from topiary._private import check
from topiary._private import installed
import pandas as pd
import subprocess, sys, os, random, string, re, warnings

def align(input_seqs,
          output_fasta=None,
          super5=False,
          silent=False,
          muscle_cmd_args=[],
          muscle_binary="muscle"):
    """
    Run muscle to align sequences.

    Parameters
    ----------
    input_seqs : str or pandas.DataFrame
        input to align (fasta file or topiary df)
    output_fasta : str or None, default=None
        output fasta file to store alignment. Optional if the input
        is a dataframe; required if the input is a fasta file.
    super5 : bool
        Use the 'super5' mode of muscle 5
    silent : bool, default=False
        whether or not to suppress all output
    muscle_cmd_args : list
        list of arguments to pass directly to muscle. Wrapper specifies :code:`-align`
        and :code:`-output` (or :code:`-in/-out` for old version of the command
        line), but leaves rest as default. Format should be something like
        :code:`['-replicates','20',...]`. Arguments are not checked by this
        function, but passed directly to muscle.
    muscle_binary : str
        location of muscle binary (default assumes 'muscle' command is in the
        :code:`$PATH`).


    Returns
    -------
    alignment : None or pandas.DataFrame
        If input_seqs is a topiary dataframe, return a copy of the dataframe
        with the aligned sequences loaded into the `alignment` column.
        Otherwise, write to `output_fasta file` and return None from the
        function.
    """

    # If output_fasta is defined, make sure it's a string
    if output_fasta and type(output_fasta) is not str:
        err = "\noutput_fasta '{output_fasta}' should be a string indicating\n"
        err += "fasta file where alignment should be written out.\n"
        raise ValueError(err)

    # If the input is a string.
    if type(input_seqs) is str:

        # Make sure the input file exists.
        if not os.path.exists(input_seqs):
            err = f"\n'{input_seqs}' file does not exist.\n\n"
            raise FileNotFoundError(err)

        # Make sure the output file is not None.
        if output_fasta is None:
            err = "\nIf aligning a fasta file, an output_fasta file must be\n"
            err += "specified.\n"
            raise ValueError(err)

        # Do the alignment.
        _run_muscle(input_seqs,output_fasta,super5,silent,muscle_cmd_args,muscle_binary)

        print(f"\nSuccess. Alignment written to '{output_fasta}'.",flush=True)

        return None

    # If the input is a dataframe
    elif type(input_seqs) is pd.DataFrame:

        # Check the dataframe
        df = check.check_topiary_dataframe(input_seqs)

        # Create temporary input file
        tmp_file_root = "".join([random.choice(string.ascii_letters) for i in range(10)])
        input_fasta = "topiary-tmp_{}_align-in.fasta".format(tmp_file_root)
        topiary.write_fasta(df,input_fasta,sort_on_taxa=False)

        # Create temporary output file name if required.
        temporary_output = False
        if output_fasta is None:
            output_fasta = "topiary-tmp_{}_align-out.fasta".format(tmp_file_root)
            temporary_output = True

        # Do the alignment
        _run_muscle(input_fasta,output_fasta,super5,silent,muscle_cmd_args,muscle_binary)

        # Read alignment back into the dataframe
        df = topiary.read_fasta_into(df,output_fasta)

        # Delete temporary input file
        try:
            os.remove(input_fasta)
        except FileNotFoundError:
            pass

        # Delete temporary output file
        if not silent:
            print("\nSuccess. Alignment written to the `alignment` column in the dataframe.",flush=True)
        if temporary_output:
            try:
                os.remove(output_fasta)
            except FileNotFoundError:
                pass
        else:
            if not silent:
                print(f"\nSuccess. Alignment written to '{output_fasta}'.",flush=True)

        return df

    else:
        err = "\ninput_seqs should either be a fasta file or a topiary dataframe\n\n"
        raise ValueError(err)


def _run_muscle(input_fasta,
                output_fasta,
                super5=False,
                silent=False,
                muscle_cmd_args=[],
                muscle_binary="muscle"):
    """
    Run muscle for sequence alignment.

    Parameters
    ----------
    input_fasta : str
        input fasta file to align
    output_fasta : str
        output fasta file to store alignment
    super5 : bool
        Use the 'super5' mode of muscle 5
    silent : bool, default=False
        whether or not to suppress all output
    muscle_cmd_args : list
        list of arguments to pass directly to muscle. Wrapper specifies :code:`-align`
        and :code:`-output` (or :code:`-in/-out` for old version of the command
        line), but leaves rest as default. Format should be something like
        :code:`['-replicates','20',...]`. Arguments are not checked by this
        function, but passed directly to muscle.
    muscle_binary : str
        location of muscle binary (default assumes 'muscle' command is in the
        :code:`$PATH`).

    Returns
    -------
    None
    """

    # Check muscle version
    binary, muscle_version = installed.check_muscle()
    if muscle_version == (-2,-2,-2):
        err = f"\nmuscle not found in the PATH ({os.environ['PATH']})\n\n"
        raise RuntimeError(err)

    if muscle_version == (-1,-1,-1):
        ret = subprocess.run([muscle_binary],capture_output=True)
        err = f"\nmuscle is in the PATH, but is crashing. The output of \n"
        err += f"`muscle` follows:\n\n {ret.stderr.decode()}\n\n"
        raise RuntimeError(err)

    if muscle_version == (0,0,0):
        w = "\nmuscle is in the PATH, but topiary could not determine the\n"
        w += "muscle version. Proceeding with assumption it is >= 5.0.\n"
        warnings.warn(w)

        muscle_version = ("5",)

    # Construct command line appropriate to muscle version being used
    if int(muscle_version[0]) >= 5:
        cmd = [muscle_binary,"-align",input_fasta,"-output",output_fasta]
        if super5:
            cmd[1] = "-super5"
    else:
        cmd = [muscle_binary,"-in",input_fasta,"-out",output_fasta]

    # Append user-specified arguments
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

        if not silent:
            # This bit of trickery makes it so the muscle output counts up to 100%
            # on a single line, rather than spewing out over 100s of lines over the
            # course of the alignment
            endl = "\n"
            if topiary._in_notebook:
                if not re.search("100.0%",line):
                    print(90*" ",end="\r")
                    endl = "\r"
            print(line.strip(),end=endl,flush=True)

    popen.stdout.close()
    return_code = popen.wait()

    # If running muscle failed...
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, cmd)
