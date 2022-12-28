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

from topiary.draw import PrettyTree

import pandas as pd
import numpy as np

import os
import inspect
import urllib
import re
import copy
import shutil

def _load_check_html(some_html):
    """
    Check validity of a block of html and return an HTMLValidator object.
    """

    parser = HTMLValidator()
    parser.feed(some_html)
    assert parser.is_valid

    return parser


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
    _load_check_html(f"{s}{e}")

    os.chdir(cwd)


def test_df_to_table():
    
    test_df = pd.DataFrame({"test":np.array([1,2,3],dtype=int),
                            "this":np.array([4.1,5.2,6.3],dtype=float)})

    # Test show_row_numbers and add_header (which interact)

    for header in [True,False]:

        out = df_to_table(test_df,
                          add_header=header,
                          show_row_numbers=True,
                          float_fmt="{}",
                          int_fmt="{}")

        parser = _load_check_html(out)

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
        parser = _load_check_html(out)

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
    parser = _load_check_html(out)
    
    assert parser.tag_dict["td"][-1][-1] == "6.30"
    assert parser.tag_dict["td"][-2][-1] == "03"

    out = df_to_table(test_df,
                      add_header=True,
                      show_row_numbers=True,
                      float_fmt="{:.3f}",
                      int_fmt="{:03d}")
    parser = _load_check_html(out)
    
    assert parser.tag_dict["td"][-1][-1] == "6.300"
    assert parser.tag_dict["td"][-2][-1] == "003"


def test_canvas_to_html():
    
    # Read as string
    tree = "(((A:1.0,B:4.0)AB:1.0,((C:1.0,D:1.0)CD:1.0,E:1.0)CDE:1.0)ABCDE:1.0,(F:1.0,G:1.0)FG)ABCDEFG;"
    pt = PrettyTree(T=tree)
    out = canvas_to_html(pt.canvas)

    _load_check_html(out)

def test_sequence_box():
    
    out = sequence_box(text="THIS IS TEST TEXT Z",
                       color="#000000",
                       prop_value=None,
                       prop_span=None,
                       palette=None)

    parser = _load_check_html(out)

    # Check some of the text inputs -- should all be in one span because single
    # color
    assert parser.tag_dict["span"][0][2] == 'THIS IS TEST TEXT Z'

    # Make sure color mapping works
    out = sequence_box(text="ABC",
                       color=["#ffffff","black"],
                       prop_value=[0,0.5,1],
                       prop_span=[0,1],
                       palette=None)

    parser = _load_check_html(out)

    assert parser.tag_dict["span"][0][2] == 'A'
    assert parser.tag_dict["span"][0][1]["style"][0] == 'color:rgba(100.0%,100.0%,100.0%,1.000)'
    assert parser.tag_dict["span"][1][2] == 'B'
    assert parser.tag_dict["span"][1][1]["style"][0] == 'color:rgba(50.0%,50.0%,50.0%,1.000)'
    assert parser.tag_dict["span"][2][2] == 'C'
    assert parser.tag_dict["span"][2][1]["style"][0] == 'color:rgba(0.0%,0.0%,0.0%,1.000)'


def test_create_card():
    
    # Defaults
    out = create_card(card_title=None,
                      card_contents=None,
                      title_tag="h6",
                      match_height=True)
    parser = _load_check_html(out)

    assert parser.tag_dict["div"][0][1]["class"][0] == "card"
    assert parser.tag_dict["div"][0][1]["class"][1] == "h-100"
    assert "h6" not in parser.tag_dict

    # match height to False
    out = create_card(card_title=None,
                      card_contents=None,
                      title_tag="h6",
                      match_height=False)
    parser = _load_check_html(out)

    assert parser.tag_dict["div"][0][1]["class"][0] == "card"
    assert len(parser.tag_dict["div"][0][1]["class"][0])
    assert "h6" not in parser.tag_dict

    # Add card title
    out = create_card(card_title="test",
                      card_contents=None,
                      title_tag="h6",
                      match_height=False)
    parser = _load_check_html(out)

    assert parser.tag_dict["div"][0][1]["class"][0] == "card"
    assert len(parser.tag_dict["div"][0][1]["class"][0])
    assert "h6" in parser.tag_dict
    assert parser.tag_dict["h6"][0][2] == "test"

    # Card title, altered title_tag
    out = create_card(card_title="test",
                      card_contents=None,
                      title_tag="h4",
                      match_height=False)
    parser = _load_check_html(out)

    assert parser.tag_dict["div"][0][1]["class"][0] == "card"
    assert len(parser.tag_dict["div"][0][1]["class"][0])
    assert "h4" in parser.tag_dict
    assert parser.tag_dict["h4"][0][2] == "test"

    # Card title, altered title_tag
    out = create_card(card_title="test",
                      card_contents="my contents",
                      title_tag="h4",
                      match_height=False)
    parser = _load_check_html(out)

    # Make sure card contents goes in
    assert parser.tag_dict["div"][0][1]["class"][0] == "card"
    assert len(parser.tag_dict["div"][0][1]["class"][0])
    assert "h4" in parser.tag_dict
    assert parser.tag_dict["h4"][0][2] == "test"

    found = False
    for div in parser.tag_dict["div"]:
        if div[2] == "my contents":
            found = True
            break
    assert found

def test_create_element():
    
    s, e = create_element("rocket")
    parser = _load_check_html(f"{s}{e}")
    
    assert len(parser.tag_dict["rocket"][0][1]) == 0
    
    s, e = create_element("rocket",{"motor":"test"})
    parser = _load_check_html(f"{s}{e}")
    
    assert np.array_equal(parser.tag_dict["rocket"][0][1]["motor"],["test"])

    s, e = create_element("rocket",{"motor":["test","this"]})
    parser = _load_check_html(f"{s}{e}")
    
    assert np.array_equal(parser.tag_dict["rocket"][0][1]["motor"],["test","this"])

def test_create_icon_row():
    
    out = create_icon_row(["stupid.txt"],["some words"])
    parser = _load_check_html(out)
    
    assert parser.tag_dict["a"][0][1]["href"][0] == "stupid.txt"
    assert np.array_equal(parser.tag_dict["a"][0][1]["title"],["some","words"])
    assert parser.tag_dict["img"][0][1]["src"][0] == ".assets/txt_icon.svg"

    out = create_icon_row(["stupid.csv"],["some words"])
    parser = _load_check_html(out)
    
    assert parser.tag_dict["a"][0][1]["href"][0] == "stupid.csv"
    assert np.array_equal(parser.tag_dict["a"][0][1]["title"],["some","words"])
    assert parser.tag_dict["img"][0][1]["src"][0] == ".assets/csv_icon.svg"

    out = create_icon_row(["stupid.pdf"],["some words"])
    parser = _load_check_html(out)
    
    assert parser.tag_dict["a"][0][1]["href"][0] == "stupid.pdf"
    assert np.array_equal(parser.tag_dict["a"][0][1]["title"],["some","words"])
    assert parser.tag_dict["img"][0][1]["src"][0] == ".assets/pdf_icon.svg"

    out = create_icon_row(["stupid.newick"],["some words"])
    parser = _load_check_html(out)
    
    assert parser.tag_dict["a"][0][1]["href"][0] == "stupid.newick"
    assert np.array_equal(parser.tag_dict["a"][0][1]["title"],["some","words"])
    assert parser.tag_dict["img"][0][1]["src"][0] == ".assets/newick_icon.svg"

    out = create_icon_row(["stupid.fasta"],["some words"])
    parser = _load_check_html(out)
    
    assert parser.tag_dict["a"][0][1]["href"][0] == "stupid.fasta"
    assert np.array_equal(parser.tag_dict["a"][0][1]["title"],["some","words"])
    assert parser.tag_dict["img"][0][1]["src"][0] == ".assets/fasta_icon.svg"

    out = create_icon_row(["stupid.doc"],["some words"])
    parser = _load_check_html(out)
    
    assert parser.tag_dict["a"][0][1]["href"][0] == "stupid.doc"
    assert np.array_equal(parser.tag_dict["a"][0][1]["title"],["some","words"])
    assert parser.tag_dict["img"][0][1]["src"][0] == ".assets/txt_icon.svg"

    out = create_icon_row(["stupid.txt","junk.csv"],["some words","more words"])
    parser = _load_check_html(out)

    assert parser.tag_dict["a"][0][1]["href"][0] == "stupid.txt"
    assert np.array_equal(parser.tag_dict["a"][0][1]["title"],["some","words"])
    assert parser.tag_dict["img"][0][1]["src"][0] == ".assets/txt_icon.svg"

    assert parser.tag_dict["a"][1][1]["href"][0] == "junk.csv"
    assert np.array_equal(parser.tag_dict["a"][1][1]["title"],["more","words"])
    assert parser.tag_dict["img"][1][1]["src"][0] == ".assets/csv_icon.svg"


def test_create_row():
    
    out = create_row(["A","B"])
    parser = _load_check_html(out)
    assert parser.tag_dict["div"][1][1]["class"][0] == "row"
    assert parser.tag_dict["div"][2][1]["class"][0] == "col"
    assert parser.tag_dict["div"][2][2] == "A"
    assert parser.tag_dict["div"][3][1]["class"][0] == "col"
    assert parser.tag_dict["div"][3][2] == "B"


def test_create_modal():
    
    out = create_modal(modal_text="modal_text",
                       modal_title="modal_title",
                       modal_label="modal_label")
    parser = _load_check_html(out)

def test_create_info_modal():
    
    out = create_info_modal(modal_text="modal_text",
                            modal_title="modal_title",
                            button_label="button_label")
    parser = _load_check_html(out)
