__description__ = \
"""
Functions for reading and writing dataframes.
"""
__author__ = "Michael J. Harms"
__date__ = "2021-04-08"

import topiary
from topiary import _arg_processors

import pandas as pd
import numpy as np

import os

def read_dataframe(input,remove_extra_index=True):
    """
    Read a topiary spreadsheet. Handles .csv, .tsv, .xlsx/.xls. If extension is
    not one of these, attempts to parse text as a spreadsheet using
    pandas.read_csv(sep=None).

    Parameters
    ----------
        input: either a pandas dataframe OR the filename to read in.
        remove_extra_index: look for the 'Unnamed: 0' column that pandas writes out
            for pd.to_csv(index=True) and, if found, drop column.

    Return
    ------
        validated topiary dataframe
    """

    # If this is a string, try to load it as a file
    if type(input) is str:

        filename = input

        ext = filename.split(".")[-1].strip().lower()

        if ext in ["xlsx","xls"]:
            df = pd.read_excel(filename)
        elif ext == "csv":
            df = pd.read_csv(filename,sep=",")
        elif ext == "tsv":
            df = pd.read_csv(filename,sep="\t")
        else:
            # Fall back -- try to guess delimiter
            df = pd.read_csv(filename,sep=None,engine="python")

    # If this is a pandas dataframe, work in a copy of it.
    elif type(input) is pd.DataFrame:
        df = input.copy()

    # Otherwise, fail
    else:
        err = f"\n\n'input' {input} not recognized. Should be the filename of\n"
        err += "spreadsheet or a pandas dataframe.\n"
        raise ValueError(err)

    # Look for extra index column that pandas writes out (in case user wrote out
    # pandas frame manually, then re-read). Looks for first column that is
    # Unnamed and has values [0,1,2,...,L]
    if remove_extra_index:
        if df.columns[0].startswith("Unnamed:"):
            possible_index = df.loc[:,df.columns[0]]
            if np.issubdtype(possible_index.dtypes,int):
                if np.array_equal(possible_index,np.arange(len(possible_index),dtype=int)):
                    df = df.drop(columns=[df.columns[0]])

    # Validate the dataframe
    df = _arg_processors.process_topiary_dataframe(df)

    return df

def write_dataframe(df,out_file,overwrite=False):
    """
    Write a dataframe to an output file. The type of file written depends on the
    extension of out_file. If .csv, write comma-separated. If .tsv, write tab-
    separated. If .xlsx, write excel. Otherwise, write as a .csv file.

    Parameters
    ----------
        df: topiary dataframe
        out_file: output file name
        overwrite: whether or not to overwrite an existing file

    Return
    ------
        None. Writes to out_file
    """

    # Validate dataframe before writing out
    df = topiary._arg_processors.process_topiary_dataframe(df)

    # Validate output file
    if type(out_file) is not str:
        err = f"\n\nout_file '{out_file}' should be a string.\n"
        raise ValueError(err)

    ext = out_file.split(".")[-1]
    if ext not in ["csv","tsv","xlsx"]:
        print("\n\nOutput file extension not recognized. Will write as csv.\n\n")
        ext = "csv"

    if os.path.exists(out_file):
        if not overwrite:
            err = f"\nout_file '{out_file}' exists. To overwrite, set overwrite to True\n\n"
            raise FileExistsError(err)
        else:
            os.remove(out_file)

    # Write out appropriate file type
    if ext == "csv":
        df.to_csv(out_file,sep=",",index=False)
    elif ext == "tsv":
        df.to_csv(out_file,sep="\t",index=False)
    else:
        df.to_excel(out_file,index=False)
