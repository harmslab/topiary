"""
Functions for reading and writing alignments to files.
"""

import topiary
from topiary import check, _private

import pandas as pd
import numpy as np

import os, re

def _validate_seq_writer(df,
                         out_file,
                         seq_column,
                         label_columns,
                         write_only_keepers,
                         empty_char,
                         clean_sequence,
                         overwrite):
    """
    Generic function for validating arguments used to validate sequences written
    out to files (such as fasta or phy).

    For argument parameter meanings, see write_fasta and write_phy.
    """

    df = topiary.check.check_topiary_dataframe(df)

    # Validate output file
    if type(out_file) is not str:
        err = f"\n\nout_file '{out_file}' should be a string.\n"
        raise ValueError(err)

    # Make sure seq column is sane
    try:
        if type(seq_column) is dict:
            raise TypeError
        df.loc[:,seq_column]
    except (KeyError,TypeError):
        err = f"seq_column '{seq_column}' not found\n."
        raise ValueError(err)

    # Make sure label columns are sane
    if hasattr(label_columns,"__iter__") and type(label_columns) is not type:
        if type(label_columns) is str:
            label_columns = [label_columns]
    else:
        err = "\nlabel_columns should be a list of column names in the\n"
        err += "dataframe.\n\n"
        raise ValueError(err)

    # Make sure the label columns are in the dataframe
    for c in label_columns:
        try:
            df[c]
        except (KeyError,TypeError):
            err = f"label_column '{c}' not recognized."
            err += "Should be a a column name in the dataframe.\n"
            raise ValueError(err)

    # Make sure label columns has uid in the first position
    try:
        label_columns.remove("uid")
    except ValueError:
        pass
    label_columns.insert(0,"uid")

    # Validate write_only_keepers argument
    try:
        bool(int(write_only_keepers))
    except (TypeError,ValueError):
        err = f"\nwrite_only_keepers '{write_only_keepers}' not recognized.\n"
        err += "Should be True or False.\n\n"
        raise ValueError(err)

    # Validate empty_char argument
    if empty_char is None:
        empty_char = []
    else:
        if type(empty_char) is not str:
            err = f"\nempty_char '{empty_char}' not recognized. This should be\n"
            err += "a string of characters\n\n"
            raise ValueError(err)

        empty_char = list(empty_char)

    # Validate clean_sequence argument
    try:
        bool(int(clean_sequence))
    except (TypeError,ValueError):
        err = f"\nclean_sequence '{clean_sequence}' not recognized. Should be\n"
        err += "True or False.\n\n"
        raise ValueError(err)

    # Validate overwrite argument
    try:
        bool(int(overwrite))
    except (TypeError,ValueError):
        err = f"\noverwrite '{overwrite}' not recognized. Should be\n"
        err += "True or False.\n\n"
        raise ValueError(err)

    # Check for file existance if overwrite = False
    if os.path.isfile(out_file):
        if not overwrite:
            err = f"\nout_file '{out_file}' already exists. To overwrite, set\n"
            err += "overwrite = True.\n\n"
            raise FileExistsError(err)

    return df, label_columns, empty_char


def write_fasta(df,out_file,seq_column="sequence",label_columns=["species","name"],
                write_only_keepers=True,empty_char="X-?",clean_sequence=False,
                overwrite=False):
    """
    Write a fasta file from a dataframe.

    Parameters
    ----------
    df: pandas.DataFrame
        data frame to write out
    out_file : str
        output file
    seq_column : str, default="sequence"
        column in data frame to use as sequence
    label_columns : list, default=["species","name"]
        list of columns to use for sequence labels
    write_only_keepers : bool, default=True
        whether or not to write only seq with keep == True
    empty_char : str or None
        string containing empty char. If the sequence is only empty char, do not
        write out. To disable check, set empty_char=None.
    clean_sequence : bool, default=False
        replace any non-aa characters with "-"
    overwrite : bool, default=False
        whether or not to overwrite an existing file

    Returns
    -------
    None
    """

    df, label_columns, empty_char = _validate_seq_writer(df,
                                                         out_file,
                                                         seq_column,
                                                         label_columns,
                                                         write_only_keepers,
                                                         empty_char,
                                                         clean_sequence,
                                                         overwrite)

    # Construct fasta output
    out = []
    for i in range(len(df)):
        row = df.iloc[i]

        if write_only_keepers:
            if not row.keep:
                continue

        # Create label for header
        h = []
        for c in label_columns:
            h.append(f"{row[c]}")
        h = "|".join(h)

        seq = row[seq_column]
        is_empty = len([s for s in seq if s not in empty_char]) == 0
        if seq == "" or seq is None or is_empty:
            continue

        # Replace non-aa characters with '-'
        if clean_sequence:
            seq = re.sub("[^ACDEFGHIKLMNPQRSTVWY-]","-",seq)

        out.append(f">{h}\n{seq}\n")

    if len(out) == 0:
        err = "\nNo sequences written out. Check your settings and the dataframe.\n\n"
        raise RuntimeError(err)

    # Write output
    f = open(out_file,"w")
    f.write("".join(out))
    f.close()

def write_phy(df,
              out_file,
              seq_column="alignment",
              write_only_keepers=True,
              empty_char="X-?",
              clean_sequence=False,
              overwrite=False):
    """
    Write a .phy file from a dataframe. Uses the uid as the sequence name. All
    sequences must have the same length.

    Parameters
    ----------
    df: pandas.DataFrame
        data frame to write out
    out_file : str
        output file
    seq_column : str, default="alignment"
        column in data frame to use as sequence
    label_columns : list, default=["species","name"]
        list of columns to use for sequence labels
    write_only_keepers : bool, default=True
        whether or not to write only seq with keep == True
    empty_char : str or None
        string containing empty char. If the sequence is only empty char, do not
        write out. To disable check, set empty_char=None.
    clean_sequence : bool, default=False
        replace any non-aa characters with "-"
    overwrite : bool, default=False
        whether or not to overwrite an existing file

    Returns
    -------
    None
    """

    label_columns = ["uid"]

    df, label_columns, empty_char = _validate_seq_writer(df,
                                                         out_file,
                                                         seq_column,
                                                         label_columns,
                                                         write_only_keepers,
                                                         empty_char,
                                                         clean_sequence,
                                                         overwrite)

    if write_only_keepers:
        num_to_write = np.sum(df.keep)
    else:
        num_to_write = len(df.keep)

    all_lengths = []
    for i in range(len(df)):
        try:
            # If this is not being kept, don't record
            if write_only_keepers and np.logical_not(df.keep.iloc[i]):
                raise TypeError

            # If sequence is length 0 (or some other wacky non-len compatible
            # input) don't record.
            l = len(df[seq_column].iloc[i])
            if l == 0:
                raise TypeError

            # Record length
            all_lengths.append(l)
        except TypeError:
            if not df.keep.iloc[i] and write_only_keepers:
                pass
            else:
                err = "\nRow does not have sequence\n\n"
                err += f"...{df.iloc[i]}\n"
                raise ValueError(err)

    all_lengths = set(all_lengths)
    if len(all_lengths) == 0:
        err = "\nno sequences found\n\n"
        raise ValueError(err)

    if len(all_lengths) != 1:
        err = "\nnot all rows have the same sequence length\n\n"
        raise ValueError(err)

    # Finally, get length of alignment
    ali_length = list(all_lengths)[0]

    # Construct phy output
    out = []
    for i in range(len(df)):
        row = df.iloc[i]

        if write_only_keepers:
            if not row.keep:
                continue

        h = row["uid"]
        seq = row[seq_column]
        is_empty = len([s for s in seq if s not in empty_char]) == 0
        if seq == "" or seq is None or is_empty:
            num_to_write -= 1
            continue

        # Replace non-aa characters with '-'
        if clean_sequence:
            seq = re.sub("[^ACDEFGHIKLMNPQRSTVWY-]","-",seq)

        out.append(f"{h}\n{seq}\n")

    if len(out) == 0:
        err = "\nNo sequences written out. Check your settings and the dataframe.\n\n"
        raise RuntimeError(err)


    out.insert(0,f"{num_to_write}  {ali_length}\n\n")

    # Write output
    f = open(out_file,"w")
    f.write("".join(out))
    f.close()

def read_fasta_into(df,fasta_file,load_into_column="alignment",unkeep_missing=True):
    """
    Load sequences from a fasta file into an existing topiary dataframe. This
    function expects the fasta file to have names formated like >uid|other stuff.
    It will match the uid in the fasta file with the uid in the topiary dataframe.
    If a uid is not in the dataframe, the function will raise an error. 

    Parameters
    ----------
    df : pandas.DataFrame
        topiary data frame
    fasta_file : str
        a fasta file with headers formatted like >uid|other stuff
    load_into_column : str, default="alignment"
        what column in the dataframe to load the sequences into
    unkeep_missing : bool, default=True
        set any sequences in the dataframe tht are not in the fasta file to
        keep=False. This allows the user to delete sequences from the alignment
        and have that reflected in the dataframe.

    Return
    ------
    pandas.DataFrame
        topiary dataframe with sequences now in load_into_column
    """

    # Create data frame and make sure it has the column in which to load
    new_df = check.check_topiary_dataframe(df)
    try:
        new_df[load_into_column]
    except KeyError:
        new_df[load_into_column] = None

    # Go through the fasta file and get sequences
    headers = []
    uids = []
    seqs = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                headers.append(line.strip())
                uids.append(line[1:].split("|")[0].strip())
                seqs.append([])
            else:
                seqs[-1].append(line.strip())

    if len(uids) != len(set(uids)):
        err = "Not all uids unique in this fasta file\n"
        raise ValueError(err)

    final_seqs = []
    for s in seqs:
        final_seqs.append("".join(s))
    uid_to_index = dict(zip(df.uid,df.index))

    # Load sequences from fasta into data frame
    loaded_seq = {}
    for i in range(len(seqs)):

        try:
            index = uid_to_index[uids[i]]
            new_df.loc[index,load_into_column] = final_seqs[i]
            loaded_seq[index] = None
        except KeyError:
            err = f"Could not map the sequence titled {headers[i]} to an index\n"
            err += "in the data frame. This function expects the sequence titles\n"
            err += "in a fasta file to have the format:\n"
            err += ">uid|other stuff here\n"
            err += f"The parsed uid ({uids[i]}) is not in the dataframe!\n\n"
            raise ValueError(err)

    # If requested, set all sequences not in alignment to Keep = False
    if unkeep_missing:
        for i in list(new_df.index):
            try:
                loaded_seq[i]
            except KeyError:
                new_df.loc[i,"keep"] = False
                new_df.loc[i,load_into_column] = pd.NA

    # Return dataframe with final sanity check to make sure uid stayed unique
    return check.check_topiary_dataframe(new_df)
