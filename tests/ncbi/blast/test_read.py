import pytest

from topiary.ncbi.blast.read import _xml_file_to_records
from topiary.ncbi.blast.read import _clean_xml
from topiary.ncbi.blast.read import check_for_cpu_limit
from topiary.ncbi.blast.read import records_to_df
from topiary.ncbi.blast.read import read_blast_xml

import os, shutil, re

def test__clean_xml(user_xml_files):

    for f in user_xml_files:
        if os.path.split(f)[-1] == "ZO1_vertebrate_nr-clustered.xml":
            test_file = f
            break

    g = open(test_file)
    lines = g.readlines()
    g.close()

    # assert this is a bad xml file
    assert re.search("CREATE_VIEW","".join(lines)) is not None

    cleaned = _clean_xml(f)
    assert re.search("CREATE_VIEW",cleaned) is None


def test__xml_file_to_records(user_xml_files):

    # Basically wrapper for biopython reader. Make sure it reads some real
    # files, including one with CREATE_VIEW nastiness
    for f in user_xml_files:

        out = _xml_file_to_records(f)
        expected_length = user_xml_files[f]["length"]
        expected_queries = user_xml_files[f]["queries"]
        assert len(out) == expected_queries

        # length is for first query
        for o in out[:1]:
            assert len(list(o.alignments)) == expected_length

def test_check_for_cpu_limit(xml):

    with pytest.raises(FileNotFoundError):
        check_for_cpu_limit("stupid")

    assert not check_for_cpu_limit(xml["good.xml"])
    assert check_for_cpu_limit(xml["cpu-limit.xml"])
    assert not check_for_cpu_limit(xml["nr_clustered.xml"])

def test_records_to_df(user_xml_files):

    for f in user_xml_files:

        xml_records = _xml_file_to_records(f)
        expected_length = user_xml_files[f]["length"]
        expected_queries = user_xml_files[f]["queries"]

        df = records_to_df(xml_records)

        # Hack, really. If there is only one query, make sure the dataframe
        # has expected length. If more than one, make sure dataframe is longer.
        # (To implement correctly -- passing full length -- I'd have to
        # refactor conftest entry. Maybe later)
        if expected_queries == 1:
            assert len(df) == expected_length
        else:
            assert len(df) > expected_length

def test_read_blast_xml(xml,tmpdir,user_xml_files):

    # Pass in a single xml file, not in a list
    xml_file = xml["good.xml"]
    df, xml_files = read_blast_xml(xml_file)
    assert isinstance(df,list)
    assert isinstance(xml_files,list)
    assert len(df) == 1
    assert len(xml_files) == 1
    assert len(df[0]) == 19
    assert xml_files[0] == xml_file

    # Pass two xml files
    df, xml_files = read_blast_xml([xml_file,xml_file])
    assert isinstance(df,list)
    assert isinstance(xml_files,list)
    assert len(df) == 2
    assert len(xml_files) == 2
    assert len(df[0]) == 19
    assert len(df[1]) == 19
    assert xml_files[0] == xml_file
    assert xml_files[1] == xml_file

    # Pass directory with an xml file
    xml_files_dir = os.path.join(tmpdir,"xml_files")
    os.mkdir(xml_files_dir)
    shutil.copy(xml_file,
                os.path.join(xml_files_dir,os.path.split(xml_file)[-1]))
    df, xml_files = read_blast_xml(xml_files_dir)

    assert len(df) == 1
    assert len(xml_files) == 1
    assert len(df[0]) == 19
    assert xml_files[0] == os.path.join(xml_files_dir,os.path.split(xml_file)[-1])

    # Make directory with no xml files. Should die.
    no_xml_files_dir = os.path.join(tmpdir,"no_xml_files")
    os.mkdir(no_xml_files_dir)
    with pytest.raises(ValueError):
        df, xml_files = read_blast_xml(no_xml_files_dir)

    # Passing a stupid xml file with broken format ... not a great test, but
    # should at least throw error of some sort.
    xml_file = xml["bad.xml"]
    with pytest.raises(ValueError):
        df, xml_files = read_blast_xml(xml_file)

    # Test read of some basic examples
    df, xml_files = read_blast_xml(user_xml_files)
    assert len(df) == len(user_xml_files)

    # Validate do_cpu_check flag by sending in something that hit a cpu limit
    # (cpu-limit.xml) and then something that did not (good.xml). If we set
    # do_cpu_check AND the file has a cpu-limit, we should see df is None.
    # Otherwise, it should not be None. 
    df, xml_files = read_blast_xml([xml["cpu-limit.xml"]],do_cpu_check=False)
    assert df is not None

    df, xml_files = read_blast_xml([xml["cpu-limit.xml"]],do_cpu_check=True)
    assert df is None

    df, xml_files = read_blast_xml([xml["good.xml"]],do_cpu_check=False)
    assert df is not None

    df, xml_files = read_blast_xml([xml["good.xml"]],do_cpu_check=True)
    assert df is not None
