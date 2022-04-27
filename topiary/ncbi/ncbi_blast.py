__author__ = "Michael J. Harms"
__date__ = "2021-04-23"
__description__ = \
"""
BLAST against an NCBI database.
"""

from . base import read_blast_xml

from Bio import Entrez
from Bio.Blast import NCBIXML, NCBIWWW
Entrez.email = "DUMMY_EMAIL@DUMMY_URL.COM"

import sys, urllib, http

def ncbi_blast(sequence,
               db="nr",
               taxid=None,
               blast_program="blastp",
               hitlist_size=50,
               e_value_cutoff=0.01,
               gapcosts=(11,1),
               num_tries_allowed=5,
               url_base="https://blast.ncbi.nlm.nih.gov/Blast.cgi",
               test_run=False,
               **kwargs):
    """
    Perform a blast query against a remote NCBI database. Takes a sequence or
    list of sequences and returns a list of topiary dataframes containing hits
    for each sequence.

    sequence: sequence as a string OR list of string sequences
    db: NCBI blast database
    taxid: taxid for limiting blast search (default to no limit)
    blast_program: NCBI blast program to use (blastp, tblastn, etc.)
    hitlist_size: download only the top hitlist_size hits
    e_value_cutoff: only take hits with e_value better than e_value_cutoff
    gapcost: BLAST gapcosts (length 2 tuple of ints)
    num_tries_allowed: try num_tries_allowed times in case of timeout
    url_base: NCBI base url
    test_run: if True, parse arguments etc. but do not actually query blast
              server. return the blast_kwargs as a dictionary.
    kwargs: extra keyword arguments are passed directly to NCBIWWW.qblast,
            overriding anything constructed above. You could, for example, pass
            entrez_query="txid9606[ORGN] or txid10090[ORGN]" to limit blast
            results to hits from human or mouse. This would take precedence over
            a taxid specified above.

    returns: dataframe (single query) or list of dataframes (multiple sequences)
             if no hits found, warns and returns None.
    """

    # -------------------------------------------------------------------------
    # Deal with sequence input

    # Make sure sequence is an iterable
    return_singleton = False
    if hasattr(sequence,"__iter__") and type(sequence) is not type:

        # If sequence is not a string, make a single fasta-formatted string
        # from contents
        out = []
        if type(sequence) is not str:
            for i, s in enumerate(sequence):
                if type(s) is not str:
                    err = "\nsequence must either be a string or list of strings, one for each sequence\n\n"
                    raise ValueError(err)

                seq = s.strip()
                if len(seq) == 0:
                    err = "\nempty sequences are not allowed\n"
                    raise ValueError(err)

                out.append(f">count{i}\n{s}\n")

            sequence = "".join(out)

        else:

            sequence = sequence.strip()
            if len(sequence) == 0:
                err = "\nempty sequences are not allowed\n"
                raise ValueError(err)

            return_singleton = True
            out.append(f">count0\n{sequence}\n")
    else:
        err = "\nsequence must either be a string or list of strings, one for each sequence\n\n"
        raise ValueError(err)

    sequence_input = "".join(out).strip()

    # -------------------------------------------------------------------------
    # Deal with db input
    if type(db) is not str:
        err = "\ndb should be a string indicating an NCBI database (e.g. 'nr')\n\n"
        raise ValueError(err)

    # -------------------------------------------------------------------------
    # Deal with taxid input

    # Convert integer input to list of one string
    if type(taxid) is int:
        taxid = [f"{taxid}"]

    taxid_is_bad = True
    taxid_out = []
    if hasattr(taxid,"__iter__") and type(taxid) is not type:

        taxid_is_bad = False

        # If taxid is a single string, convert it to a list of one string
        if type(taxid) is str:
            taxid = [taxid]

        # Go through list of taxids passed in
        for t in taxid:
            if type(t) in [int,str]:
                taxid_out.append(f"txid{t}[ORGN]")
            else:
                taxid_is_bad = True
                break

    else:
        if taxid is None:
            taxid_out = []
            taxid_is_bad = False

    if taxid_is_bad:
        err = "\ntaxid should be either a single ncbi taxon id (e.g. 9606 for\n"
        err += "human) or list of ncbi taxon ids. Individual ids can either be\n"
        err += "int or str. These are passed to NCBI without validation.\n\n"
        raise ValueError(err)

    # -------------------------------------------------------------------------
    # Construct qblast keywords

    # keyword arguments to pass to qblast
    qblast_kwargs = {"program":blast_program,
                     "database":db,
                     "sequence":sequence_input,
                     "hitlist_size":f"{hitlist_size}",
                     "expect":f"{e_value_cutoff}",
                     "gapcosts":f"{gapcosts[0]} {gapcosts[1]}",
                     "url_base":url_base}

    # Construct taxid entrez_query
    if len(taxid_out) > 0:
        qblast_kwargs["entrez_query"] = " or ".join(taxid_out)

    # Capture the rest of kwargs, overwriting anything automatically made
    for k in kwargs:
        qblast_kwargs[k] = kwargs[k]

    # If a test run, do not actually do blast, but instead return blast kwargs.
    if test_run:
        return qblast_kwargs

    # -------------------------------------------------------------------------
    # Run blast query

    # While we haven't run out of tries...
    tries = 0
    while tries < num_tries_allowed:

        print(f"Running NCBI blast against {db}. Attempt {tries+1} of {num_tries_allowed}.")
        sys.stdout.flush()

        # Try to read output.
        try:

            result = NCBIWWW.qblast(**qblast_kwargs)
            out = NCBIXML.parse(result)

        # If some kind of http error or timeout, set out to None
        except (urllib.error.URLError,urllib.error.HTTPError,http.client.IncompleteRead):
            out = None

        # If out is None, try again. If not, break out of loop--success!
        if out is None:
            tries += 1
            continue
        else:
            break

    # We didn'g et result even after num_tries_allowed tries. Throw
    # an error.
    if out is None:
        err = "\nProblem accessing with NCBI server\n\n"
        raise RuntimeError(err)

    print("Success.")
    sys.stdout.flush()

    out_df = read_blast_xml(out)

    # If no hits were found, return None
    if len(out_df) == 0:
        w = "\nNo hits returned. This can happen if there is a silent error\n"
        w += "on a remote server. Try changing your query and blasting again.\n\n"
        warnings.warn(w)
        return None

    # If user passed in single sequence (not list) return a single dataframe
    # instead of a list of dataframes
    if return_singleton:
        out_df = out_df[0]

    return out_df
