#!/usr/bin/env python3
"""
Command line interface to topiary.io.alignments.read_fasta_into
"""

import topiary
from topiary.io.alignments import read_fasta_into
from topiary._private.wrap import wrap_function

import sys

def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    # Set arg types for args with None as default
    optional_arg_types = {}

    extra_args = [("out",{"type":str})]

    description = \
    """
    Load sequences from a fasta file into an existing topiary dataframe. This
    function expects the fasta file to have names formated like >uid|other stuff.
    It will match the uid in the fasta file with the uid in the topiary dataframe.
    If a uid is not in the dataframe, the function will raise an error.

    Parameters
    ----------
    df : str
        file with topiary dataframe
    fasta : str
        a fasta file with headers formatted like >uid|other stuff
    out : str
        output file
    load_into_column : str, default="alignment"
        what column in the dataframe to load the sequences into
    unkeep_missing : bool, default=True
        set any sequences in the dataframe tht are not in the fasta file to
        keep=False. This allows the user to delete sequences from the alignment
        and have that reflected in the dataframe.

    Return
    ------
    pandas.DataFrame
        topiary dataframe with sequences now in load_into_column
    """

    # Wrap and run function
    ret, args = wrap_function(read_fasta_into,
                              argv=argv,
                              optional_arg_types=optional_arg_types,
                              extra_args=extra_args,
                              description=description)

    out_file = args.__dict__["out"]
    topiary.write_dataframe(ret,out_file)

if __name__ == "__main__":
    main()
