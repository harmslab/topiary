__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
Input/output functions for topiary.
"""

from . import ncbi
from . import util
from . import _private
from . import opentree

import pandas as pd
import numpy as np

# Modules for blasting, etc.
from Bio import SeqIO

from tqdm.auto import tqdm

import re, string, random, io, warnings, copy




def read_dataframe(input,remove_extra_index=True):
    """
    Read a spreadsheet. Handles .csv, .tsv, .xlsx/.xls. If extension is not one
    of these, attempts to parse text as a spreadsheet using
    pandas.read_csv(sep=None).

    input: either a pandas dataframe OR the filename to read in.
    remove_extra_index: look for the 'Unnamed: 0' column that pandas writes out
        for pd.to_csv(index=True) and, if found, drop column.
    """

    # If this is a string, try to load it as a file
    if type(input) is str:

        filename = input

        ext = filename.split(".")[-1]

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
    df = util.check_topiary_dataframe(df)

    return df

def write_dataframe(df,out_file):
    """
    Write a dataframe to an output file. The type of file written depends on the
    extension of out_file. If .csv, write comma-separated. If .tsv, write tab-
    separated. If .xlsx, write excel. Otherwise, write as a .csv file.

    df: topiary dataframe
    out_file: output file name
    """

    if type(out_file) is not str:
        err = f"\n\nout_file '{out_file}' should be a string.\n"
        raise ValueError(err)

    ext = out_file.split(".")[-1]
    if ext not in ["csv","tsv","xlsx"]:
        w = "\n\nOutput file extension not recognized. Will write as csv.\n\n"
        warnings.warn(w)
        ext = "csv"

    # Write out appropriate file type
    if ext == "csv":
        df.to_csv(out_file,sep=",",index=False)
    elif ext == "tsv":
        df.to_csv(out_file,sep="\t",index=False)
    else:
        df.to_excel(out_file,index=False)



def ncbi_blast_xml_to_df(xml_files):
    """
    Take a list of blast xml files and load in all sequences as a single
    topiary data frame. Parse meta data in an intelligent way, download
    sequences via entrez, and find unique taxonomic identifiers on the open
    tree of life.

    xml_files: blast xml files to load. if a string, treat as a single xml file.
               if a list, treat as a list of xml files.

    returns a pandas data frame
    """

    # If only one xml file is specified, convert it to a list (of one) xml
    # file
    if type(xml_files) is str:
        xml_files = [xml_files]

    # List to hold all hits and accession numbers to download
    all_hits = []
    to_download = []

    # For each xml file
    for i, xml in enumerate(xml_files):

        # Read xml
        tmp_df = ncbi.read_blast_xml(xml)

        # Go through and parse meta data
        for j in range(len(tmp_df)):

            accession = tmp_df.loc[j,"accession"]
            title = tmp_df.loc[j,"title"]
            evalue = tmp_df.loc[j,"e_value"]

            start = tmp_df.loc[j,"subject_start"]
            end = tmp_df.loc[j,"subject_end"]
            hit_info = ncbi.parse_ncbi_line(title,
                                            accession=accession)
            if hit_info is None:
                continue

            all_hits.append((xml,start,end,evalue,hit_info))
            to_download.append(hit_info["accession"])


    # Reduce to a unique set of accessions
    if len(set(to_download)) != len(to_download):
        unique_hits = []
        unique_download = []
        for i, d in enumerate(to_download):
            if d in unique_download:
                continue

            unique_hits.append(all_hits[i])
            unique_download.append(d)

        all_hits = unique_hits
        to_download = unique_download

    # Download sequences from entrez
    all_output = ncbi.entrez_download(to_download)

    # Capture sequences from the downloaded data
    captured = []
    for record in SeqIO.parse(io.StringIO(all_output), "fasta"):
        seq_id = str(record.id)
        sequence = str(record.seq)
        captured.append((seq_id,sequence))

    # Create a dictionary with appropriate keys to load into dataframe.
    out = util.create_pipeline_dict()

    # Add blast-specific columns
    blast_key_list = ["accession","xml","length","evalue","start","end",
                      "structure","low_quality","precursor","predicted",
                      "isoform","hypothetical","partial","raw_line"]
    for k in blast_key_list:
        out[k] = []

    # Go through every hit
    for i in range(len(all_hits)):

        # Get information from previous few rounds
        seq = captured[i][1]
        accession = captured[i][0]

        xml = all_hits[i][0]
        start = all_hits[i][1]
        end = all_hits[i][2]
        evalue = all_hits[i][3]
        hit_info = all_hits[i][4]

        # Get hit_info, if k is in key_list
        for k in hit_info.keys():
            try:
                out[k].append(hit_info[k])
            except KeyError:
                pass

        # Overwrite accession from hit_info
        out["accession"][-1] = accession

        # Load info from blast itself
        out["xml"].append(xml)
        out["sequence"].append(seq)
        out["length"].append(len(seq))
        out["evalue"].append(evalue)
        out["start"].append(start)
        out["end"].append(end)

        out["uid"].append(_private.generate_uid())

        out["keep"].append(True)

    df = pd.DataFrame(out)
    df = util.check_topiary_dataframe(df)

    return df


def write_fasta(df,out_file,seq_column="sequence",seq_name="pretty",
                write_only_keepers=True,empty_char="X-?",clean_sequence=False):
    """
    df: data frame to write out
    out_file: output file
    seq_column: column in data frame to use as sequence
    seq_name: column in data frame to use as >NAME.  If "pretty",
              write out a pretty names.
    write_only_keepers: whether or not to write only seq with keep = True
    empty_char: empty char. if the sequence is only empty char, do not write
                out.
    clean_sequence: replace any non-aa characters with "-"
    """

    # Make sure seq name is sane
    try:
        df[seq_name]
        take_pretty = False
    except KeyError:
        if seq_name == "pretty":
            take_pretty = True
        else:
            err = f"seq_name '{seq_name}' not recognized."
            err += "Should be a column name or 'pretty'\n"
            raise ValueError(err)

    # Make sure seq column is sane
    try:
        df[seq_column]
    except KeyError:
        err = f"seq_column '{seq_column}' not found\n."
        raise ValueError(err)

    # Construct fasta output
    out = []
    for i in range(len(df)):
        row = df.iloc[i]

        if write_only_keepers:
            if not row.keep:
                continue

        if take_pretty:
            h = _private.to_pretty(row)
        else:
            h = row[seq_name]

        seq = row[seq_column]
        is_empty = len([s for s in seq if s not in list(empty_char)]) == 0
        if seq == "" or seq is None or is_empty:
            continue

        # Replace non-aa characters with '-'
        if clean_sequence:
            seq = re.sub("[^ACDEFGHIKLMNPQRSTVWYZ-]","-",seq)

        out.append(f">{h}\n{seq}\n")

    # Write output
    f = open(out_file,"w")
    f.write("".join(out))
    f.close()


def write_phy(df,out_file,seq_column="sequence",
              write_only_keepers=True,
              empty_char="X-?",
              clean_sequence=False):
    """
    Write out a .phy file using uid as keys.

    df: data frame to write out
    out_file: output file
    seq_column: column in data frame to use as sequence
    write_only_keepers: whether or not to write only seq with keep = True
    empty_char: empty char. if the sequence is only empty char, do not write
                out.
    clean_sequence: replace any non-aa characters with "-"
    """

    # Make sure seq column is sane
    try:
        df[seq_column]
    except KeyError:
        err = f"seq_column '{seq_column}' not found\n."
        raise ValueError(err)

    if write_only_keepers:
        num_to_write = np.sum(df.keep)
    else:
        num_to_write = len(df.keep)

    all_lengths = []
    for i in range(len(df)):
        try:
            l = len(df[seq_column].iloc[i])
            if l == 0:
                raise TypeError
            all_lengths.append(l)
        except TypeError:
            if not df.keep.iloc[i] and write_only_keepers:
                pass
            else:
                err = "\n\nRow does not have sequence\n\n"
                err += f"...{df.iloc[i]}\n"
                raise ValueError(err)

    all_lengths = set(all_lengths)
    if len(all_lengths) == 0:
        err = "\n\nno sequences found\n"
        raise ValueError(err)

    if len(all_lengths) != 1:
        err = "\n\nnot all rows have the same alignment length\n\n"
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
        is_empty = len([s for s in seq if s not in list(empty_char)]) == 0
        if seq == "" or seq is None or is_empty:
            num_to_write -= 1
            continue

        # Replace non-aa characters with '-'
        if clean_sequence:
            seq = re.sub("[^ACDEFGHIKLMNPQRSTVWYZ-]","-",seq)

        out.append(f"{h}\n{seq}\n")

    out.insert(0,f"{num_to_write}  {ali_length}\n\n")

    # Write output
    f = open(out_file,"w")
    f.write("".join(out))
    f.close()

def load_fasta(df,fasta_file,load_into_column="alignment",empty_char="X-?",unkeep_missing=True):
    """
    unkeep_missing: set sequences not loading into keep = False
    empty_char: empty char. if the sequence is only empty char, set keep = False
    """

    # Create new data frame and make sure it has the column in which to load
    new_df = df.copy()
    try:
        new_df[load_into_column]
    except KeyError:
        new_df[load_into_column] = None

    # Figure out how pretty and uid calls map to df index
    pretty_to_index, uid_to_index = _private.get_index_maps(new_df)

    # Go through the fasta file and get sequences
    header = []
    seqs = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                if len(seqs) > 0:
                    seqs[-1] = "".join(seqs[-1])
                header.append(line.strip()[1:])
                seqs.append([])
            else:
                seqs[-1].extend(list(line.strip()))
    seqs[-1] = "".join(seqs[-1])

    # Load sequences from fasta into data frame
    loaded_seq = {}
    for i in range(len(header)):

        # Figure out the index to modify
        try:
            index = pretty_to_index[header[i]]
        except KeyError:
            try:
                index = uid_to_index[header[i]]
            except KeyError:
                err = f"could not map {header[i]} to data frame\n"
                raise ValueError(err)

        # Actually modify data frame
        new_df.loc[index,load_into_column] = seqs[i]

        # Record the sequence was loaded if it's not all junk (like -?X)
        if len([s for s in seqs[i] if s not in list(empty_char)]) > 0:
            loaded_seq[index] = None

    # If requested, set all sequences not in alignment to Keep = False
    if unkeep_missing:
        for i in list(new_df.index):
            try:
                loaded_seq[i]
            except KeyError:
                new_df.loc[i,"keep"] = False

    return new_df
