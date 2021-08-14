
import pytest
import topiary

# ---------------------------------------------------------------------------- #
# Test __init__
# ---------------------------------------------------------------------------- #

def test_create_pipeline_dict():

    pipe_dict = topiary.util.create_pipeline_dict()

    key_list = ["accession","protein","species","xml","sequence",
                "length","evalue","start","end",
                "structure","low_quality","precursor","predicted","isoform",
                "hypothetical","partial","raw_line","uid","keep"]

    assert set(pipe_dict.keys()) == set(key_list)

def test_grab_line_meta_data(ncbi_lines,ncbi_lines_parsed):

    for i, line in enumerate(ncbi_lines):
        line_dict = topiary.util.grab_line_meta_data(line)

        out = []
        for k in ncbi_lines_parsed[i]:
            try:
                out.append(line_dict[k] == ncbi_lines_parsed[i][k])
            except KeyError:
                pass

        assert sum(out) == len(ncbi_lines_parsed[i])

#def pretty_to_uid(df,to_convert,out_file=None,overwrite=False):
#def uid_to_pretty(df,to_convert,out_file=None,overwrite=False):
#def load_tree(tree,fmt=None):
