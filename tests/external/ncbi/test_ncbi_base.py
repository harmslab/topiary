
import pytest

import topiary
from topiary.external.ncbi.base import _standard_blast_args_checker as _sbac

import copy

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

def test__standard_blast_args_checker():

    default_kwargs = {"sequence":"test",
                      "hitlist_size":50,
                      "e_value_cutoff":0.01,
                      "gapcosts":(11,1)}

    kwargs = copy.deepcopy(default_kwargs)

    kwargs["sequence"] = "test"
    sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)
    assert len(sequence_list) == 1
    assert sequence_list[0] == "test"
    assert return_singleton == True

    kwargs["sequence"] = ["test"]
    sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)
    assert len(sequence_list) == 1
    assert sequence_list[0] == "test"
    assert return_singleton == False

    kwargs["sequence"] = ["test","this"]
    sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)
    assert len(sequence_list) == 2
    assert sequence_list[0] == "test"
    assert sequence_list[1] == "this"
    assert return_singleton == False

    # Send in bad sequences
    bad_sequences = [1,False,1.0,-1,None,str,""]
    for b in bad_sequences:
        print("passing bad sequence:",b)

        kwargs["sequence"] = b
        with pytest.raises(ValueError):
            sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)

        kwargs["sequence"] = [b]
        with pytest.raises(ValueError):
            sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)

    # -------------------------------------------------------------------------
    # histlist_size

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)
    assert hitlist_size == 50

    kwargs["hitlist_size"] = 100
    sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)
    assert hitlist_size == 100

    bad_int = [0,False,[],-1,None,str,"",{},int]
    for b in bad_int:
        print("passing bad hitlist_size:",b)
        kwargs["hitlist_size"] = b
        with pytest.raises(ValueError):
            sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)



    # -------------------------------------------------------------------------
    # evalue_cutoff

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)
    assert e_value_cutoff == 0.01

    kwargs["e_value_cutoff"] = 0.001
    sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)
    assert e_value_cutoff == 0.001

    bad_float = [-1,[],None,str,"",{},int]
    for b in bad_float:

        kwargs["e_value_cutoff"] = b
        print("passing bad e_value_cutoff:",b)
        with pytest.raises(ValueError):
            sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)

    # -------------------------------------------------------------------------
    # gapcosts

    kwargs = copy.deepcopy(default_kwargs)

    sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)
    assert gapcosts == (11,1)

    kwargs["gapcosts"] = [10,1]
    sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)
    assert gapcosts == (10,1)

    bad_gap = [-1,[],[1,2,3],[-1,1],[1,-1],["test","this"],None,str,"",{},int]
    for b in bad_gap:
        kwargs["gapcosts"] = b
        print("passing bad gapcosts:",b)
        with pytest.raises(ValueError):
            sequence_list, hitlist_size, e_value_cutoff, gapcosts, return_singleton = _sbac(**kwargs)
