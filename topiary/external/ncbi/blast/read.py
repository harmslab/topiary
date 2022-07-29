"""
Read BLAST xml output.
"""

import numpy as np
import pandas as pd

from Bio.Blast import NCBIXML

import io, os, glob


def _xml_file_to_records(xml_file):
    """
    Convert the contents of an xml file into a list of blast records. There will
    be one record per query. This function cleans up the "CREATE_VIEW" mangling
    that ncbi sometimes injects into xml.

    Parameters
    ----------
    xml_file : str
        xml file to load

    Returns
    -------
    blast_records : list
        list of blast records, one per query in xml file
    """

    # Read xml file, stripping mangled blank lines and CREATE_VIEW
    # that NCBI sometimes injects into otherwise valid xml
    lines = []
    with open(xml_file) as f:
        for line in f:
            if line.strip() not in ["","CREATE_VIEW"]:
                lines.append(line)
    file_contents = "".join(lines)

    # Now use Biopython to parse blast records
    blast_records = []
    with io.StringIO(file_contents) as f:
        for rec in NCBIXML.parse(f):
            blast_records.append(rec)

    return blast_records

def records_to_df(blast_records):
    """
    Read list of blast records and return as a single pandas data frame.

    Parameters
    ----------
    blast_records : list
        list of biopython.Blast records

    Returns
    -------
    results : list
        pandas dataframe with all blast hits
    """

    out_df = []
    for record in blast_records:

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

def read_blast_xml(xml_input):
    """
    Load blast xml file(s) and convert to pandas dataframe(s).

    Parameters
    ----------
    xml_input : str or list
        blast xml files to load. This can have a few formats:

        1) single xml file (str)
        2) list of xml files
        3) directory (str). topiary will grab all .xml in that directory.

    Returns
    -------
    all_df : list
        pandas dataframes constructed from the blast xml files.
    xml_files : list
        list of xml files parsed
    """

    xml_files = []
    if isinstance(xml_input,str):

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

        # Make sure it is an iterable (type check catches edge case where user
        # sends in the list datatype)
        if not hasattr(xml_input,"__iter__") or isinstance(xml_input,type):
            err = f"\nCould not parse xml_input '{xml_input}'. Should be an xml\n"
            err += "file, list of xml files, or directory containing .xml files.\n\n"
            raise ValueError(err)

        for x in xml_input:
            if os.path.isfile(x):
                xml_files.append(x)
            else:
                err = f"\nxml file '{x}' not found\n"
                raise ValueError(err)

    # Actually parse xml files
    all_df = []
    for x in xml_files:
        records = _xml_file_to_records(x)
        all_df.append(records_to_df(records))

    return all_df, xml_files
