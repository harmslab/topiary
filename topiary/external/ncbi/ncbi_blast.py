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

import numpy as np
import pandas as pd

import sys, urllib, http, copy

def ncbi_blast(sequence,
               db="nr",
               taxid=None,
               blast_program="blastp",
               hitlist_size=50,
               e_value_cutoff=0.01,
               gapcosts=(11,1),
               num_tries_allowed=5,
               max_query_length=80000,
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
    max_query_length: maximum string length accepted by the server. if the
                      query is too long, this function will break it into
                      multiple requests, each sent to ncbi.
    url_base: NCBI base url
    test_run: if True, parse arguments etc. but do not actually query blast
              server. return the blast_kwargs as a dictionary.
    kwargs: extra keyword arguments are passed directly to Bio.Blast.NCBIWWW.qblast,
            overriding anything constructed above. You could, for example, pass
            entrez_query="txid9606[ORGN] or txid10090[ORGN]" to limit blast
            results to hits from human or mouse. This would take precedence over
            any taxid specified above.

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

    db_is_bad = True
    if type(db) is str:
        if len(db) > 1:
            db_is_bad = False

    if db_is_bad:
        err = "\ndb should be a string indicating an NCBI database (e.g. 'nr')\n\n"
        raise ValueError(err)

    # -------------------------------------------------------------------------
    # Deal with taxid input

    taxid_out = []
    taxid_is_bad = True

    # Is it int-like? Convert integer input to list of one string. We're
    # stringent on this -- no floats -- because an int cast will silently
    # round down. If someone used a taxid like 960.6 (extra .) it would be
    # interpreted as 960, which is very much the wrong taxid
    if np.issubdtype(type(taxid),np.integer):
        taxid = [f"{taxid}"]

    # This wacky line sees if the taxid is iterable, but not a type *class*.
    # Catches weird edge case where user passes in str or int as a taxid
    if hasattr(taxid,"__iter__") and type(taxid) is not type:

        taxid_is_bad = False

        # If taxid is a single string, convert it to a list of one string
        if type(taxid) is str:
            taxid = [taxid]

        # Go through list of taxids and put in correct format
        for t in taxid:
            if type(t) is str:
                taxid_out.append(f"txid{t}[ORGN]")
            elif np.issubdtype(type(t),np.integer):
                taxid_out.append(f"txid{t}[ORGN]")
            else:
                taxid_is_bad = True
                break

    # If taxid was None to begin with, ignore it
    else:
        if taxid is None:
            taxid_out = []
            taxid_is_bad = False

    if taxid_is_bad:
        err = "\ntaxid should be either a single ncbi taxon id (e.g. 9606 for\n"
        err += "human) or list of ncbi taxon ids. Individual ids can either be\n"
        err += "int or str. These are passed to NCBI without further validation.\n\n"
        raise ValueError(err)

    # -------------------------------------------------------------------------
    # Deal with blast_program

    blast_program_is_bad = True
    if type(blast_program) is str:
        if len(blast_program) > 1:
            blast_program_is_bad = False

    if blast_program_is_bad:
        err = "\nblast_program should be a string identifier for the blast\n"
        err += "program to use (i.e. blastp, tblastn, etc.)\n"
        raise ValueError(err)

    # -------------------------------------------------------------------------
    # Deal with hitlist_size

    try:
        hitlist_size = int(hitlist_size)
        if hitlist_size < 1:
            raise ValueError
    except (ValueError,TypeError):
        err = "\nhitlist_size should be an integer greater than one.\n"
        raise ValueError(err)

    # -------------------------------------------------------------------------
    # Deal with e_value_cutoff

    try:
        e_value_cutoff = float(e_value_cutoff)
        if e_value_cutoff < 0:
            raise ValueError
    except (ValueError,TypeError):
        err = "\ne_value_cutoff should be a positive float value\n"
        raise ValueError(err)

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

    # -------------------------------------------------------------------------
    # Deal with num_tries_allowed

    try:
        num_tries_allowed = int(num_tries_allowed)
        if num_tries_allowed < 1:
            raise ValueError
    except (ValueError,TypeError):
        err = "\num_tries_allowed should be an integer greater than zero.\n"
        raise ValueError(err)

    # -------------------------------------------------------------------------
    # Deal with max_query_length

    try:
        max_query_length = int(max_query_length)
        if max_query_length < 1:
            raise ValueError
    except (ValueError,TypeError):
        err = "\nmax_query_length should be an integer greater than one.\n"
        raise ValueError(err)

    # -------------------------------------------------------------------------
    # Deal with url_base

    url_base_is_bad = True
    if type(url_base) is str:
        if len(url_base) > 1:
            url_base_is_bad = False

    if url_base_is_bad:
        err = "\nurl_base should be a string holding a url pointing to the\n"
        err += "blast server you wish to use.\n"
        raise ValueError(err)


    # -------------------------------------------------------------------------
    # Construct qblast keywords

    # keyword arguments to pass to qblast *besides* sequence
    qblast_kwargs = {"program":blast_program,
                     "database":db,
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

    # -------------------------------------------------------------------------
    # Deal with sequence kwargs

    # check to see if the sequence input is too long
    if len(sequence_input)//(max_query_length + 1) > 0:

        # Split sequences up
        seqs = sequence_input.split("\n>")

        # List will queries of useful lengths
        tmp_out = [[]]
        current_split_size = 0
        for s in seqs:
            if len(s) + current_split_size < max_query_length:
                tmp_out[-1].append(s)
                current_split_size += len(s)
            else:

                # If you have a super, super long sequence...
                if len(s) > max_query_length:
                    err = "\nat least one sequence is, by itself, greater than\n"
                    err += f"the maximum server query length '{max_query_length}'.\n"
                    err += "To fix, either trim this sequence, remove it from the\n"
                    err += "query list, increase the max_query_length (not\n"
                    err += "recommended), or use a local blast database rather\n"
                    err += "than the NCBI database.\n"
                    raise ValueError(err)

                # Rip off trailing "\n>"
                tmp_out[-1][-1] = tmp_out[-1][-1][:-2]

                # Create new
                tmp_out.append([])
                tmp_out[-1].append(s)
                current_split_size = len(s)

        # Recompile sequence chunks into fasta queries.
        split_seqs = []
        for t in tmp_out:
            split_seqs.append("\n>".join(t))

            # Every sequence but the first will need a new > added to start
            if not split_seqs[-1].startswith(">"):
                split_seqs[-1] = f">{split_seqs[-1]}"

        sequence_input = split_seqs

        print(f"Broke blast query into {len(sequence_input)} queries\n")
        sys.stdout.flush()

    else:
        sequence_input = [sequence_input]

    # Now make a list of queries to make
    queries_to_make = []
    for s in sequence_input:
        queries_to_make.append(copy.deepcopy(qblast_kwargs))
        queries_to_make[-1]["sequence"] = s

    # If a test run, do not actually do blast, but instead return blast kwargs.
    if test_run:
        return queries_to_make

    # -------------------------------------------------------------------------
    # Run blast query

    print(f"Running NCBI blast against {db}.")
    out_df = []
    for query_counter, this_query in enumerate(queries_to_make):

        print(f"    Query {query_counter + 1} of {len(queries_to_make)}.")
        sys.stdout.flush()

        # While we haven't run out of tries...
        tries = 0
        while tries < num_tries_allowed:

            print(f"        Attempt {tries+1} of {num_tries_allowed}.")
            sys.stdout.flush()

            # Try to read output.
            try:

                result = NCBIWWW.qblast(**this_query)
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

        # We didn't get result even after num_tries_allowed tries. Throw
        # an error.
        if out is None:
            err = "\nProblem accessing with NCBI server. We got no output after\n"
            err += f"{num_tries_allowed} attempts. This is likely due to a server\n"
            err += "timeout. Try again at a different time.\n"
            raise RuntimeError(err)

        print(f"    Successfully completed query {query_counter +1} of {len(queries_to_make)}")
        sys.stdout.flush()

        # NCBI returns a single blast record, with all query sequences in a single
        # big dataframe.
        ncbi_df = read_blast_xml(out)

        # Break big dataframe into a list of dataframes, one for each query sequence
        queries = np.unique(ncbi_df["query"])
        query_order = [(q[5:],q) for q in queries]
        query_order.sort()
        for q in query_order:

            this_df = ncbi_df.loc[ncbi_df["query"] == q[1],:]

            # No hits, return empty dataframe
            if len(this_df) == 1 and pd.isna(this_df["accession"].iloc[0]):
                out_df.append(pd.DataFrame())

            # Record hits
            else:
                out_df.append(this_df)

    print("Done blasting.")
    sys.stdout.flush()

    # If no hits were found, return None
    if len(out_df) == 0:
        w = "\nNo hits returned. This can happen if there is a silent error\n"
        w += "on a remote server. Try changing your query and blasting again.\n"
        print(w)
        return None

    # If user passed in single sequence (not list) return a single dataframe
    # instead of a list of dataframes
    if return_singleton:
        out_df = out_df[0]

    return out_df
