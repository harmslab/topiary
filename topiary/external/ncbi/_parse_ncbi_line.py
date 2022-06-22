"""
Functions for interacting with NCBI databases and file types.
"""

from topiary import util, check

import pandas as pd
import numpy as np

# Modules for blasting, etc.
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.Blast import NCBIXML, NCBIWWW
import Bio.Blast.Applications as apps
from Bio import pairwise2

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

    Parameters
    ----------
        line: line from an NCBI record of some sort

    Return
    ------
        dictionary with true or false for different patterns in line
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


def parse_ncbi_line(line,accession=None):
    """
    Parse an ncbi line of the sort seen in the BLAST title field or on each
    line of a fasta file.

    Parameters
    ----------
        line: line from an NCBI record
        accession: extract entry from line that matches acccession.  Ignores
                   version (e.g. "1" in XXXXXXXXX.1).  If None, parse first
                   entry on line.

    Return
    ------
        dictionary with following keys:
            raw_line: unprocessed line (input)
            line: processed line (remove multiple titles)
            name: protein name
            structure: whether or not this is a structure (bool)
            low_quality: whether or not this is low quality (bool)
            predicted: whether or not this is predicted (bool)
            precursor: whether or not this is a precursor (bool)
            isoform: whether or not this is an isoform (bool)
            hypothetical: whether or not this is hypothetical (bool)
            partial: whether or not this is a partial sequence (bool)
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
