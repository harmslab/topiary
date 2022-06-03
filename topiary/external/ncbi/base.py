__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
functions interacting with ncbi databases and file types.
"""

from topiary import util

import pandas as pd
import numpy as np

# Modules for blasting, etc.
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.Blast import NCBIXML, NCBIWWW
import Bio.Blast.Applications as apps
from Bio import pairwise2

Entrez.email = "DUMMY_EMAIL@DUMMY_URL.COM"

from tqdm.auto import tqdm

import re, sys, os, string, random, pickle, io, urllib, http, subprocess
import multiprocessing as mp

def _grab_line_meta_data(line):
    """
    Look for meta data we care about from the line.  This includes:

        + structure
        + low_quality
        + precursor
        + predicted
        + isoform
        + hypothetical
        + partial

    Return a dictionary with each of those keyed to a bool for whether
    or not the line has this.
    """

    meta_patterns = {"structure":"crystal structure",
                     "low_quality":"low.quality",
                     "predicted":"predicted",
                     "precursor":"precursor",
                     "isoform":"isoform",
                     "hypothetical":"hypothetical",
                     "partial":"partial"}

    out = dict([(p,None) for p in meta_patterns])
    for m in meta_patterns:
        mp = re.compile(meta_patterns[m],re.IGNORECASE)
        out[m] = bool(mp.search(line))

    return out

def read_blast_xml(blast_record):
    """
    Read BLAST XML format and return as a pandas data frame.

    blast_record: either biopython Blast record or string pointing to blast xml file

    returns: list of dataframes, one for each query in blast_record
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

def parse_ncbi_line(line,accession=None):
    """
    Parse an ncbi line of the sort seen in the BLAST title field or on each
    line of a fasta file.

    accession: extract entry from line that matches acccession.  Ignores version
               (e.g. "1" in XXXXXXXXX.1).  If None, parse first entry on line.

    Returns a dictionary with following keys:
        raw_line -> unprocessed line (input)
        line -> processed line (remove multiple titles)
        name -> protein name
        structure -> whether or not this is a structure (bool)
        low_quality -> whether or not this is low quality (bool)
        predicted -> whether or not this is predicted (bool)
        precursor -> whether or not this is a precursor (bool)
        isoform -> whether or not this is an isoform (bool)
        hypothetical -> whether or not this is hypothetical (bool)
        partial -> whether or not this is a partial sequence (bool)
    """

    out = {"raw_line":line}

    # Clean up leading and trailing spaces
    line = line.strip()

    # Clean up leading ">" if it's there (difference between fasta and blast
    # lines)
    if line.startswith(">"):
        line = line[1:]

    # Particularly in BLAST hits, there are often multiple entries separated
    # by ">" (basically, pasted together fasta entries). Only take one.  If
    # accession is specified, take that accesion.  If no accession is
    # specified, take the first one.

    # Split on ">"
    line_splitter = re.compile(">")
    entries_on_line = [s.strip() for s in line_splitter.split(line)]

    # Extract accession from entry (assumes xxx|acccession stuff).  Ignores
    # trailing accession version number XXXXXXXX.1 -> XXXXXXXX
    accession_dict = {}
    for e in entries_on_line:
        try:
            if e.startswith("pdb"):
                k = [v.split()[0] for v in e.split("|")[1:3]]
                k = "{}_{}".format(*k)
            else:
                k = e.split("|")[-2].split(".")[0].strip()
            accession_dict[k] = e
        except IndexError:
            print(f"Could not parse line {e}. Skipping.")
            return None

    # If accession is specified, grab entry corresponding to that from the
    # line.  If not, grab first entry.
    if not accession:
        accession = list(accession_dict.keys())[0]

    accession = accession.split(".")[0].strip()
    try:
        line = accession_dict[accession]
    except KeyError:
        warn = f"accession '{accession}' not found in line '{line}'\n"
        warn += f"Using accession {accession} anyway and returning entire line."
        print(warn)

    out["accession"] = accession

    # Record the processed line
    out["line"] = line

    # Grab meta stuff (predicted, isoform, etc.)
    meta = _grab_line_meta_data(line)
    for m in meta:
        out[m] = meta[m]

    # Look for species name (thing within [ xxx ])
    species = None
    species_pattern = re.compile("\[.*?]")

    sm = species_pattern.search(line)
    if sm:
        out["species"] = sm.group(0)[1:-1]
    else:
        out["species"] = None

    # Clean up any double spaces introduced into the line at this point
    line = re.sub("  "," ",line)

    # Protein name (takes between '| XXXXXXXXX [' ).
    out["name"] = line.split("|")[-1].split("[")[0].strip()

    return out
