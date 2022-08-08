
import pytest

import topiary
import copy

from topiary.external.ncbi._parse_ncbi_line import parse_ncbi_line, _grab_line_meta_data

def test__grab_line_meta_data(ncbi_lines):

    # This function is implicitly and rigorously tested by test_parse_ncbi_line.
    # I'd have to make separate test data manually, so I'll rely on
    # test_parse_ncbi_line.
    return True


def test_parse_ncbi_line(ncbi_lines):

    input_lines = ncbi_lines[0]
    ncbi_lines_parsed = ncbi_lines[1]

    for i, line in enumerate(input_lines):
        line_dict = parse_ncbi_line(line)

        out = []
        for k in ncbi_lines_parsed[i]:

            print(line,i,k)
            assert line_dict[k] == ncbi_lines_parsed[i][k]

    print("Test ability to grab other entry off the line rather than first entry")
    line_dict = parse_ncbi_line(input_lines[0],accession="AIC51010")

    expected = copy.deepcopy(ncbi_lines_parsed[0])
    expected["line"] = "gb|AIC51010.1| LY96, partial [synthetic construct]"
    expected["accession"] = "AIC51010"
    expected["name"] = "LY96, partial"
    expected["species"] = "synthetic construct"
    expected["precursor"] = False
    expected["partial"] = True

    for k in line_dict:
        print(k)
        assert line_dict[k] == expected[k]
