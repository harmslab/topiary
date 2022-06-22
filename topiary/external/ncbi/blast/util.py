"""
Shared functions for dealing with NCBI blast.
"""

from topiary import util, check

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



def read_blast_xml(blast_record):
    """
    Read BLAST XML format and return as a pandas data frame.

    Parameters
    ----------
    blast_record : biopython.Blast record or str
        either biopython Blast record or string pointing to blast xml file

    Returns
    -------
    results : list
        list of dataframes, one for each query in blast_record
    """

    # If argument is a string, treat as xml file
    if type(blast_record) is str:

        f = open(blast_record, 'r')
        p = NCBIXML.parse(f)

        blast_record = []
        for r in p:
            blast_record.append(r)
        f.close()

    out_df = []
    for record in blast_record:

        # Prepare DataFrame fields.
        data = {'accession': [],
                'hit_def': [],
                'hit_id': [],
                'title': [],
                'length': [],
                'e_value': [],
                'sequence': [],
                'subject_start': [],
                'subject_end':[],
                'query_start':[],
                'query_end':[],
                'query':[]}

        # Get alignments from blast result.
        for i, s in enumerate(record.alignments):
            data['accession'].append(s.accession)
            data['hit_def'].append(s.hit_def)
            data['hit_id'].append(s.hit_id)
            data['title'].append(s.title)
            data['length'].append(s.length)
            data['e_value'].append(s.hsps[0].expect)
            data['sequence'].append(s.hsps[0].sbjct)
            data['subject_start'].append(s.hsps[0].sbjct_start)
            data['subject_end'].append(s.hsps[0].sbjct_end)
            data['query_start'].append(s.hsps[0].query_start)
            data['query_end'].append(s.hsps[0].query_end)
            data['query'].append(record.query)

        if len(record.alignments) == 0:
            data['accession'].append(pd.NA)
            data['hit_def'].append(pd.NA)
            data['hit_id'].append(pd.NA)
            data['title'].append(pd.NA)
            data['length'].append(pd.NA)
            data['e_value'].append(pd.NA)
            data['sequence'].append(pd.NA)
            data['subject_start'].append(pd.NA)
            data['subject_end'].append(pd.NA)
            data['query_start'].append(pd.NA)
            data['query_end'].append(pd.NA)
            data['query'].append(record.query)

        # Port to DataFrame.
        out_df.append(pd.DataFrame(data))

    out_df = pd.concat(out_df)

    return out_df
