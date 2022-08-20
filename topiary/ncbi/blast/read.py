"""
Read BLAST xml output.
"""

import numpy as np
import pandas as pd

from Bio.Blast import NCBIXML

import io, os, glob, re
import xml.etree.ElementTree as ET

def _clean_xml(xml_file):
    """
    This function cleans up the "CREATE_VIEW" mangling that ncbi sometimes
    injects into xml.

    Parameters
    ----------
    xml_file : str
        xml file to load

    Returns
    -------
    file_contents : str
        xml contents as a single string
    """

    # Read xml file, stripping mangled blank lines and CREATE_VIEW
    # that NCBI sometimes injects into otherwise valid xml
    lines = []
    with open(xml_file) as f:
        for line in f:
            if line.strip() not in ["","CREATE_VIEW"]:
                lines.append(line)

    file_contents = "".join(lines)

    return file_contents

def _xml_file_to_records(xml_file):
    """
    Convert the contents of an xml file into a list of blast records. There will
    be one record per query.

    Parameters
    ----------
    xml_file : str
        xml file to load

    Returns
    -------
    blast_records : list
        list of blast records, one per query in xml file
    """

    file_contents = _clean_xml(xml_file)

    # Now use Biopython to parse blast records
    blast_records = []
    with io.StringIO(file_contents) as f:
        for rec in NCBIXML.parse(f):
            blast_records.append(rec)

    return blast_records

def check_for_cpu_limit(xml_file):
    """
    Check to see if an ncbi server rejected the request because it hit a CPU
    limit.

    Parameters
    ----------
    xml_file : str
        xml file that came off server

    Returns
    -------
    result : bool
        True if the cpu limit was hit, False otherwise.
    """

    xml_file = str(xml_file)
    if not os.path.isfile(xml_file):
        err = f"\nxml_file '{xml_file}' does not exist.\n\n"
        raise FileNotFoundError(err)

    # This text is passed out as an <Iteration_message>
    fail_pattern = re.compile("CPU usage limit was exceeded")

    # Load file contents, making sure xml is not mangled by NCBI server
    file_contents = _clean_xml(xml_file)

    root = ET.fromstring(file_contents)
    for msg in root.iter("Iteration_message"):
        if fail_pattern.search(msg.text):
            return True

    return False

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
                'bits': [],
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
            data['bits'].append(s.hsps[0].bits)
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
            data['bits'].append(pd.NA)
            data['sequence'].append(pd.NA)
            data['subject_start'].append(pd.NA)
            data['subject_end'].append(pd.NA)
            data['query_start'].append(pd.NA)
            data['query_end'].append(pd.NA)
            data['query'].append(record.query)

        # Port to DataFrame.
        out_df.append(pd.DataFrame(data))

    out_df = pd.concat(out_df,ignore_index=True)

    return out_df

def read_blast_xml(xml_input,do_cpu_check=False):
    """
    Load blast xml file(s) and convert to pandas dataframe(s).

    Parameters
    ----------
    xml_input : str or list
        blast xml files to load. This can have a few formats:

        1) single xml file (str)
        2) list of xml files
        3) directory (str). topiary will grab all .xml in that directory.

    do_cpu_check : bool, default=False
        check files to see if they indicate cpu limit exceeded. if True and
        this is seen, return None, xml_files.

    Returns
    -------
    all_df : list or None
        pandas dataframes constructed from the blast xml files. Will be None
        if cpu_limit_check == True and the cpu limit was observed.
    xml_files : list
        list of xml files parsed.
    """

    xml_files = []
    if issubclass(type(xml_input),str):

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

    # If do_cpu_check requested, check for failed files before parsing. If we
    # indeed have a cpu failure, return None and list of xml files
    if do_cpu_check:
        for x in xml_files:
            if check_for_cpu_limit(x):
                return None, xml_files

    # Actually parse xml files
    all_df = []
    for x in xml_files:
        records = _xml_file_to_records(x)
        all_df.append(records_to_df(records))

    return all_df, xml_files
