
__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
utility functions for the topiary package
"""

from . import _private

import pandas as pd
import numpy as np

import re, sys, os, string, random, pickle, io, urllib, http

def create_pipeline_dict():
    """
    Create a dictionary with keys for column names and lists for values.
    This can be populated by a loop through sequences, followed by
    pandas.DataFrame(this_dictionary) to create a pandas data frame of
    the sort expected by the functions used in this module.
    """

    # Column names for output dictionary
    key_list = ["accession","protein","species","xml","sequence",
                "length","evalue","start","end",
                "structure","low_quality","precursor","predicted","isoform",
                "hypothetical","partial","raw_line","uid","keep"]

    return dict([(k,[]) for k in key_list])

def grab_line_meta_data(line):
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

def pretty_to_uid(some_file,out_file,df):
    """
    Use contents of data frame to convert from pretty names to uid within
    some text file.

    some_file: text file to edit
    out_file: output file name
    df: dataframe with pretty name data and uid
    """
    _private._convert_file(some_file,out_file,df,uid_to_pretty=False)

def uid_to_pretty(some_file,out_file,df):
    """
    Use contents of data frame to convert from uid to pretty names
    some text file.

    some_file: text file to edit
    out_file: output file name
    df: dataframe with pretty name data and uid
    """

    _private._convert_file(some_file,out_file,df,uid_to_pretty=True)
