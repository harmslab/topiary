import pytest

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
