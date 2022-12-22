import pytest

from conftest import HTMLValidator

from topiary.reports.elements import create_output_directory
from topiary.reports.elements import create_main_html
from topiary.reports.elements import df_to_table
from topiary.reports.elements import canvas_to_html
from topiary.reports.elements import sequence_box
from topiary.reports.elements import create_card
from topiary.reports.elements import create_element
from topiary.reports.elements import create_icon_row
from topiary.reports.elements import create_row
from topiary.reports.elements import create_modal
from topiary.reports.elements import create_info_modal

import pandas as pd

import os
import inspect
import urllib
import re
import copy
import shutil



def test_create_output_directory(tmpdir):

    cwd = os.getcwd()
    os.chdir(tmpdir)

    create_output_directory("test")
    assert os.path.isdir("test")
    assert os.path.isdir(os.path.join("test",".assets"))

    with pytest.raises(FileExistsError):
        create_output_directory("test")

    create_output_directory("test",overwrite=True)
    assert os.path.isdir("test")
    assert os.path.isdir(os.path.join("test",".assets"))

    os.chdir(cwd)

def test_create_main_html(tmpdir):
    
    cwd = os.getcwd()
    os.chdir(tmpdir)

    def _check_default(kwarg):

        # Make sure no file is there
        try:
            os.remove("default_file")
        except FileNotFoundError:
            pass

        # Get file and make sure it can still be downloaded
        sig = inspect.signature(create_main_html)
        default_file = sig.parameters[kwarg].default
        urllib.request.urlretrieve(default_file,"default_file")
        assert os.path.isfile("default_file")

    # Make sure the cdn's coded in are correct
    _check_default("bs_css")
    _check_default("bs_js")

    # Make sure bits are assembled (roughly) correctly
    s, e = create_main_html(description="TEST_AAA",
                            title="TEST_BBB",
                            custom_css="TEST_CCC",
                            bs_css="TEST_DDD",
                            bs_css_key="TEST_EEE",
                            bs_js="TEST_FFF",
                            bs_js_key="TEST_GGG")

    assert len(s.split("TEST_AAA")) == 2
    assert len(e.split("TEST_AAA")) == 1

    assert len(s.split("TEST_BBB")) == 2
    assert len(e.split("TEST_BBB")) == 1

    assert len(s.split("TEST_CCC")) == 2
    assert len(e.split("TEST_CCC")) == 1

    assert len(s.split("TEST_DDD")) == 2
    assert len(e.split("TEST_DDD")) == 1

    assert len(s.split("TEST_EEE")) == 2
    assert len(e.split("TEST_EEE")) == 1

    assert len(s.split("TEST_FFF")) == 1
    assert len(e.split("TEST_FFF")) == 2

    assert len(s.split("TEST_GGG")) == 1
    assert len(e.split("TEST_GGG")) == 2

    # This check makes sure the elements are nested correctly
    parser = HTMLValidator()
    parser.feed(f"{s}{e}")
    assert parser.is_valid

    os.chdir(cwd)


def test_df_to_table():
    
    test_df = pd.DataFrame({"test":[1,2,3],"this":[4.1,5.2,6.3]})

    # Test show_row_numbers and add_header (which interact)

    for header in [True,False]:

        out = df_to_table(test_df,
                          add_header=header,
                          show_row_numbers=True,
                          float_fmt="{}",
                          int_fmt="{}")

        parser = HTMLValidator()
        parser.feed(out)
        assert parser.is_valid

        # Look for <thead> -- only there if header is present
        if header:
            parser.tag_dict["thead"]
        else:
            with pytest.raises(KeyError):
                parser.tag_dict["thead"]

        # Look for a tag like this <th>#</th>. Will be there if header, but 
        # not if no header.
        found = False
        for t in parser.tag_dict["th"]:
            if t[2] == "#":
                found = True
                break        
        if header:
            assert found
        else:
            assert not found


        out = df_to_table(test_df,
                          add_header=header,
                          show_row_numbers=False,
                          float_fmt="{}",
                          int_fmt="{}")
        parser = HTMLValidator()
        parser.feed(out)
        assert parser.is_valid

        # Look for <thead> -- only there if header is present
        if header:
            parser.tag_dict["thead"]
        else:
            with pytest.raises(KeyError):
                parser.tag_dict["thead"]

        # Look for a tag like this <th>#</th>. Will not be there if header is 
        # True; no <th> at all if header is False
        if header:        
            for t in parser.tag_dict["th"]:
                if t[2] == "#":
                    found = True
                    break      
        else:
            assert not "th" in parser.tag_dict

        # Make sure we have correct number of rows and elements
        if header:
            assert len(parser.tag_dict["tr"]) == 4
        else:
            assert len(parser.tag_dict["tr"]) == 3
        
        if header:
            assert len(parser.tag_dict["td"]) == 6


    # Validate ability to pass float and integer format
    out = df_to_table(test_df,
                      add_header=True,
                      show_row_numbers=True,
                      float_fmt="{:.2f}",
                      int_fmt="{:02d}")
    parser = HTMLValidator()
    parser.feed(out)
    assert parser.is_valid
    
    assert parser.tag_dict["td"][-1][-1] == "6.30"
    assert parser.tag_dict["td"][-2][-1] == "03"

    out = df_to_table(test_df,
                      add_header=True,
                      show_row_numbers=True,
                      float_fmt="{:.3f}",
                      int_fmt="{:03d}")
    parser = HTMLValidator()
    parser.feed(out)
    assert parser.is_valid
    
    assert parser.tag_dict["td"][-1][-1] == "6.300"
    assert parser.tag_dict["td"][-2][-1] == "003"



def test_canvas_to_html():
    pass

def test_sequence_box():
    pass

def test_create_card():
    pass

def test_create_element():
    pass

def test_create_icon_row():
    pass

def test_create_row():
    pass

def test_create_modal():
    pass

def test_create_info_modal():
    pass