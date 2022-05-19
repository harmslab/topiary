
import topiary

def test__grab_line_meta_data(ncbi_lines,ncbi_lines_parsed):

    for i, line in enumerate(ncbi_lines):
        line_dict = topiary.external.ncbi.base._grab_line_meta_data(line)

        out = []
        for k in ncbi_lines_parsed[i]:
            try:
                out.append(line_dict[k] == ncbi_lines_parsed[i][k])
            except KeyError:
                pass

        assert sum(out) == len(ncbi_lines_parsed[i])
