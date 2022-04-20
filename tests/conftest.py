import pytest
import pandas as pd
import os, glob

@pytest.fixture(scope="module")
def ncbi_lines():
    """
    A list of generic ncbi lines.
    """

    return ["sp|Q9Y6Y9|LY96_HUMAN Lymphocyte antigen 96 OS=Homo sapiens OX=9606 GN=LY96 PE=1 SV=2"]

@pytest.fixture(scope="module")
def ncbi_lines_parsed():
    """
    A list of parse results for those ncbi_lines.
    """

    parsed_lines =  [{"structure":False,
                      "low_quality":False,
                      "predicted":False,
                      "precursor":False,
                      "isoform":False,
                      "hypothetical":False,
                      "partial":False}]

    return parsed_lines

@pytest.fixture(scope="module")
def dataframe_good_files():
    """
    A list of files that should be properly parsable by topiary.
    """

    dir = os.path.dirname(os.path.realpath(__file__))
    files = glob.glob(os.path.join(dir,"data","good","*"))

    return files



@pytest.fixture(scope="module")
def test_dataframes():
    """
    A dictionary holding dataframes of varying sorts of badness to test parser.
    """

    dir = os.path.dirname(os.path.realpath(__file__))
    search_string = os.path.join(dir,"data","*.csv")

    df_dict = {}
    for g in glob.glob(search_string):
        key = ".".join(os.path.split(g)[1].split(".")[:-1])
        df_dict[key] = pd.read_csv(g)

    return df_dict
