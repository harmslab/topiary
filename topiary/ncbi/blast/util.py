"""
Shared functions for dealing with NCBI blast.
"""

from topiary import util
from topiary._private import check

import pandas as pd
import numpy as np

# Modules for blasting, etc.
from Bio.Blast import NCBIXML

def _standard_blast_args_checker(sequence,
                                 hitlist_size,
                                 e_value_cutoff,
                                 gapcosts):
    """
    Parse and validate input used by all blast queries.

    Parameters
    ----------
        sequence: sequence as a string OR list of string sequences
        hitlist_size: download only the top hitlist_size hits
        e_value_cutoff: only take hits with e_value better than e_value_cutoff
        gapcosts: BLAST gapcosts (length 2 tuple of ints)

    Return
    ------
        list of sequence, hitlist_size, e_value_cutoff, gapcosts (tuple), and
        whether or not to return final output as a single sequence or list of
        sequences.
    """

    # Parse the sequence input. Should either be a list of strings or a single
    # string sequence. Create sequence_list which is a list of sequences (even
    # if only one)
    return_singleton = False
    if hasattr(sequence,"__iter__") and type(sequence) is not type:

        # If sequence is not a string, make a single fasta-formatted string
        # from contents
        sequence_list = []
        if type(sequence) is not str:
            for i, s in enumerate(sequence):
                if type(s) is not str:
                    err = "\nsequence must either be a string or list of strings, one for each sequence\n\n"
                    raise ValueError(err)

                seq = s.strip()
                if len(seq) == 0:
                    err = "\nempty sequences are not allowed\n"
                    raise ValueError(err)

                sequence_list.append(seq)

        else:

            sequence = sequence.strip()
            if len(sequence) == 0:
                err = "\nempty sequences are not allowed\n"
                raise ValueError(err)

            return_singleton = True
            sequence_list.append(sequence)
    else:
        err = "\nsequence must either be a string or list of strings, one for each sequence\n\n"
        raise ValueError(err)

    # -------------------------------------------------------------------------
    # Deal with hitlist_size, e_value_cutoff

    hitlist_size = check.check_int(hitlist_size,
                                               "hitlist_size",
                                               minimum_allowed=1)
    e_value_cutoff = check.check_float(e_value_cutoff,
                                                   "e_value_cutoff",
                                                   minimum_allowed=0,
                                                   minimum_inclusive=False)

    # -------------------------------------------------------------------------
    # Deal with gapcosts

    try:
        g0 = int(gapcosts[0])
        g1 = int(gapcosts[1])
        if g0 < 1 or g1 < 1:
            raise ValueError
        if len(gapcosts) > 2:
            raise ValueError

        gapcosts = (g0,g1)
    except (IndexError,ValueError,TypeError,KeyError):
        err = "\ngapcosts should be a list of two integers greater than one\n"
        raise ValueError(err)

    return sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton
