
import pandas as pd
import numpy as np

import re, sys, os, string, random, pickle, io, urllib, http

def _to_pretty(row):
    """
    Given a pandas Series, create pretty output.
    """

    try:
        pretty = f"{row.paralog}|{row.species}|{row.accession}"
    except AttributeError:
        err = "\n\nrow does not have all required attributes:"
        err += " (paralog, species, accession)\n"
        raise ValueError(err)

    return pretty

def convert(df,to_convert,out_file=None,overwrite=False,
            uid_to_pretty=False):
    """
    Private function wrapped by uid_to_pretty and pretty_to_uid that converts
    a file or string between uid and pretty.

    df: dataframe with pretty name data and uid
    to_convert: content to edit. if this is a file, read in. If not, treat as
                a text string to edit.
    out_file: output file name. If specified, write to a file.
    overwrite: if writing an output file, whether or not to overwrite.
    uid_to_pretty: if True, uid->pretty; if False, pretty->uid

    returns converted string
    """

    # If the file specifies an input file, read it in and convert that file
    if os.path.isfile(to_convert):
        f = open(to_convert,'r')
        some_string = f.read()
        f.close()
    else:
        some_string = to_convert

    # Convert the string
    for i in range(len(df)):

        # Grab uid and pretty
        row = df.iloc[i]
        uid = row.uid
        pretty = _to_pretty(row)
        escaped_pretty = re.escape(pretty)

        # Set search and convert_to
        if uid_to_pretty:
            search = re.compile("{}".format(uid))
            convert_to = pretty
        else:
            search = re.compile(escaped_pretty)
            convert_to = uid

        # Do substitutions
        some_string = search.sub(convert_to,some_string)

    # If an output file is specified
    if out_file is not None:

        if os.path.isfile(out_file):
            if not overwrite:
                err = f"file {out_file} already exists.\n"
                raise FileExistsError(err)

        f = open(out_file,'w')
        f.write(some_string)
        f.close()

    return some_string

def get_index_maps(df):
    """
    Create dictionaries that map pretty_to_uid and uid_to_pretty.

    df: data frame with columns to create pretty content and uid
    """

    pretty_to_index = {}
    uid_to_index = {}

    # Figure out pretty_to_index and uid_to_index mapping
    for i in range(len(df)):

        row = df.iloc[i,:]
        pretty = _to_pretty(row)
        uid = row.uid

        try:
            pretty_to_index[pretty]
            err = f"pretty key '{pretty}' duplicated in data frame.\n"
            raise ValueError(err)
        except KeyError:
            pretty_to_index[pretty] = i

        try:
            uid_to_index[uid]
            err = f"uid '{uid}' duplicated in data frame.\n"
            raise ValueError(err)
        except KeyError:
            uid_to_index[uid] = i

    return pretty_to_index, uid_to_index
