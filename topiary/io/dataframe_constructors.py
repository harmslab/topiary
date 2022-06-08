__description__ = \
"""
Functions for building topiary dataframes.
"""
__author__ = "Michael J. Harms"
__date__ = "2021-04-08"

import topiary
from topiary import _arg_processors, _private, ncbi

import pandas as pd
import numpy as np

from Bio import SeqIO

import os, glob, io

def ncbi_blast_xml_to_df(xml_input):
    """
    Take a list of blast xml files and load in all sequences as a single
    topiary data frame. Parse meta data in an intelligent way, download
    sequences via entrez, and find unique taxonomic identifiers on the open
    tree of life.

    Parameters
    ----------
        xml_input: blast xml files to load. This can have a few formats:
            1. single xml file
            2. list of xml files
            3. directory (grabs all files matching .xml in that directory)

    Return
    ------
        a topiary data frame
    """

    if type(xml_input) is str:

        # Looks like a file; treat as one
        if os.path.isfile(xml_input):
            xml_files = [xml_input]
        else:
            if os.path.isdir(xml_input):
                xml_files = glob.glob(os.path.join(xml_input,"*.xml"))
                xml_files.sort()
            else:
                xml_files = []

        if len(xml_files) == 0:
            err = f"\nCould not parse xml_input. Tried to read xml_input\n"
            err += f"'{xml_input}' as a file, then as a directory with .xml\n"
            err += "files. No xml files found. xml_input should be an xml file,\n"
            err += "list of xml files, or directory containing .xml files.\n\n"
            raise ValueError(err)

    else:
        xml_files = []
        if hasattr(xml_input,"__iter__") and type(xml_input) is not type:
            for x in xml_input:
                if os.path.isfile(x):
                    xml_files.append(x)
                else:
                    err = f"\nxml file '{x}' not found\n"
                    raise ValueError(err)
        else:
            err = f"\nCould not parse xml_input '{xml_input}'. Should be an xml\n"
            err += "file, list of xml files, or directory containing .xml files.\n\n"
            raise ValueError(err)


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
    out = topiary._private.create_pipeline_dict()

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
        out["start"].append(start-1)
        out["end"].append(end)

        out["uid"].append(_private.generate_uid())

        out["keep"].append(True)

    df = pd.DataFrame(out)
    df = _arg_processors.process_topiary_dataframe(df)

    return df


def fasta_to_df(fasta_input):

    # TRY NCBI, TRY UNIPROT,
    # Look for OS= OR

    # Try to parse with uniprot
    # Try to parse with ncbi

    # UNIPROT
    #>db|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier \
    #   [GN=GeneName ]PE=ProteinExistence SV=SequenceVersion
    # cols = line[1:].split("|")
    # db = cols[0]
    # accession = col1[1]
    # entry = col[1]
    # pattern = re.compile("..=")
    # for field in pattern.split(entry):
    #     print(field)

    pass
