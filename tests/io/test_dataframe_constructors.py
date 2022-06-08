import pytest

from topiary import io
import numpy as np
import pandas as pd

import warnings, os, shutil, re

def test_ncbi_blast_xml_to_df(xml,tmpdir):

    # Pass in a single xml file, not in a list
    xml_file = xml["good.xml"]
    df = io.ncbi_blast_xml_to_df(xml_file)
    assert len(df) == 19

    # Pass two xml files (indetical, so should end up with single 19-length df)
    df = io.ncbi_blast_xml_to_df([xml_file,xml_file])
    assert len(df) == 19

    # Pass directory with an xml file
    xml_files_dir = os.path.join(tmpdir,"xml_files")
    os.mkdir(xml_files_dir)
    shutil.copy(xml_file,
                os.path.join(xml_files_dir,os.path.split(xml_file)[-1]))
    df = io.ncbi_blast_xml_to_df(xml_files_dir)
    assert len(df) == 19

    # Passing a stupid xml file with broken format ... not a great test, but
    # should at least throw error of some sort.
    xml_file = xml["bad.xml"]
    with pytest.raises(ValueError):
        df = io.ncbi_blast_xml_to_df(xml_file)
