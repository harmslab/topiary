
import pandas as pd
import numpy as np

import re, sys, os, string, random, pickle, io, urllib, http

def _to_pretty(row):
    """
    Given a pandas Series, create pretty output.
    """

    try:
        pretty = f"{row.rev_call}|{row.species}|{row.accession}"
    except AttributeError:
        err = "\n\nrow does not have all required attributes:"
        err += " (rev_call, species, accession)\n"
        raise ValueError(err)

    return pretty

def _convert_file(some_file,out_file,df,uid_to_pretty=False):
    """
    Convert uid to pretty or vice versa within a file.
    Private.  Should call pretty_to_uid or uid_to_pretty wrapper functions.

    some_file: file to work on
    out_file: output
    df: dataframe with uid and contents for _to_pretty mapping
    uid_to_pretty: if True, uid->pretty; if False, pretty->uid
    """

    # Load file
    f = open(some_file)
    file_string = f.read()
    f.close()

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
        file_string = search.sub(convert_to,file_string)

    # Write output file
    f = open(out_file,'w')
    f.write(file_string)
    f.close()

def _get_index_maps(df):
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
