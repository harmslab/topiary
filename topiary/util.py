
__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
utility functions for the topiary package
"""

from . import _private

import ete3
from ete3 import Tree
import dendropy as dp

import pandas as pd
import numpy as np

import re, sys, os, string, random, pickle, io, urllib, http
import warnings

def create_pipeline_dict():
    """
    Create a dictionary with keys for column names and lists for values.
    This can be populated by a loop through sequences, followed by
    pandas.DataFrame(this_dictionary) to create a pandas data frame of
    the sort expected by the functions used in this module.
    """

    # Column names for output dictionary
    key_list = ["accession","protein","species","xml","sequence",
                "length","evalue","start","end",
                "structure","low_quality","precursor","predicted","isoform",
                "hypothetical","partial","raw_line","uid","keep"]

    return dict([(k,[]) for k in key_list])

def grab_line_meta_data(line):
    """
    Look for meta data we care about from the line.  This includes:

        + structure
        + low_quality
        + precursor
        + predicted
        + isoform
        + hypothetical
        + partial

    Return a dictionary with each of those keyed to a bool for whether
    or not the line has this.
    """

    meta_patterns = {"structure":"crystal structure",
                     "low_quality":"low.quality",
                     "predicted":"predicted",
                     "precursor":"precursor",
                     "isoform":"isoform",
                     "hypothetical":"hypothetical",
                     "partial":"partial"}

    out = dict([(p,None) for p in meta_patterns])
    for m in meta_patterns:
        mp = re.compile(meta_patterns[m],re.IGNORECASE)
        out[m] = bool(mp.search(line))

    return out

def pretty_to_uid(df,to_convert,out_file=None,overwrite=False):
    """
    Use contents of data frame to convert from pretty names to uid within
    some text file.

    df: dataframe with pretty name data and uid
    to_convert: content to edit. if this is a file, read in. If not, treat as
                a text string to edit.
    out_file: output file name. If specified, write to a file.
    overwrite: if writing an output file, whether or not to overwrite.

    returns converted string
    """

    return _private.convert(df,to_convert,out_file,overwrite,uid_to_pretty=False)


def uid_to_pretty(df,to_convert,out_file=None,overwrite=False):
    """
    Use contents of data frame to convert from uid to pretty names
    some text file.

    df: dataframe with pretty name data and uid
    to_convert: content to edit. if this is a file, read in. If not, treat as
                a text string to edit.
    out_file: output file name. If specified, write to a file.
    overwrite: if writing an output file, whether or not to overwrite.

    returns converted string
    """

    return _private.convert(df,to_convert,out_file,overwrite,uid_to_pretty=True)


def load_tree(tree,fmt=None):
    """
    Load a tree into an ete3 tree data structure.

    tree: some sort of tree. can be an ete3.Tree (returns self), a dendropy
          Tree (converts to newick and drops root), a newick file or a newick
          string.
    fmt: format for reading tree from newick.  0-9 or 100. See ete3 documentation
         for how these are read (http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees).
         As of ETE3.1.1, these numbers mean:


         |        ======  ==============================================
         |        FORMAT  DESCRIPTION
         |        ======  ==============================================
         |        0        flexible with support values
         |        1        flexible with internal node names
         |        2        all branches + leaf names + internal supports
         |        3        all branches + all names
         |        4        leaf branches + leaf names
         |        5        internal and leaf branches + leaf names
         |        6        internal branches + leaf names
         |        7        leaf branches + all names
         |        8        all names
         |        9        leaf names
         |        100      topology only
         |        ======  ==============================================

         if fmt is None, try to parse without a format descriptor, then these
         formats in numerical order.

    Returns an ete3 tree object.
    """

    # Already an ete3 tree.
    if type(tree) is ete3.TreeNode:
        return tree

    # Convert dendropy tree into newick (drop root)
    if type(tree) is dp.Tree:
        tree = tree.as_string(schema="newick",suppress_rooting=True)

    # If we get here, we need to convert. If fmt is not specified, try to parse
    # without a format string.
    if fmt is None:


        try:
            t = Tree(tree)
        except ete3.parser.newick.NewickError:

            # Try all possible formats now, in succession
            w = "\n\nCould not parse tree without format string. Going to try different\n"
            w += "formats. Please check output carefully.\n\n"
            warnings.warn(w)

            formats = list(range(10))
            formats.append(100)

            t = None
            for f in formats:
                try:
                    t = Tree(tree,format=f)
                    w = f"\n\nSuccessfully parsed tree with format style {f}.\n"
                    w += "Please see ete3 documentation for details:\n\n"
                    w += "http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees\n\n"
                    warnings.warn(w)
                    break

                except ete3.parser.newick.NewickError:
                    continue

            if t is None:
                err = "\n\nCould not parse tree!\n\n"
                raise ValueError(err)

    else:
        # Try a conversion with the specified format
        t = Tree(tree,format=fmt)

    return t
