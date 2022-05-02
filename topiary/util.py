
__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
utility functions for the topiary package
"""

from . import _private
from . import opentree

import ete3
from ete3 import Tree
import dendropy as dp

import pandas as pd
import numpy as np

import re, sys, os, string, random, pickle, io, urllib, http, copy

def create_pipeline_dict():
    """
    Create a dictionary with keys for column names and lists for values.
    This can be populated by a loop through sequences, followed by
    pandas.DataFrame(this_dictionary) to create a pandas data frame of
    the sort expected by the functions used in this module.
    """

    # Column names for output dictionary
    key_list = ["name","species","sequence","uid","keep"]

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
            print(w)

            formats = list(range(10))
            formats.append(100)

            t = None
            for f in formats:
                try:
                    t = Tree(tree,format=f)
                    w = f"\n\nSuccessfully parsed tree with format style {f}.\n"
                    w += "Please see ete3 documentation for details:\n\n"
                    w += "http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees\n\n"
                    print(w)
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

def check_topiary_dataframe(df):
    """
    Check to make sure topiary dataframe is sane. Edits dataframe, returning a
    copy.

    + Checks for required columns 'species', 'name', and 'sequence'
    + Deletes duplicated rows.
    + Deletes empty rows.

    + 'keep':
        + Must be boolean.
        + If present but not boolean, try to coerce it into boolean.
            + if str, map regex [1yt] --> True and regex[0nf] --> False
            + if int, map !0 --> True and 0 --> False
            + if float, map !~0 --> True and ~0 --> False
        + If keep is not present, add and set all to True.
    + 'uid':
        + Must be unique str of length 10 consisting only of [a-z][A-Z]
        + If present:
            + If missing values, fill them in
            + If values do not meet rules, replace them with uid that do
            + If values are not unique, update them to be unique
        + If not present, create them.
    + 'ott': If present, makes sure it has the format ottINTEGER.

    + Enforces following column order:
        + nickname (if present)
        + keep
        + species
        + name
        + sequence
        + all other columns, in order they came in
    """

    # Make sure type is right
    if type(df) is not pd.DataFrame:
        err = "\ndf must be a pandas dataframe.\n\n"
        raise ValueError(err)

    # Make sure the dataframe has all required columns
    for c in _private.required_columns:

        try:
            df[c]
        except KeyError:
            err = f"\n\nA topiary dataframe must have a {c} column.\n\n"
            raise ValueError(err)

    # Drop duplicates
    df = df.drop_duplicates(ignore_index=True)

    # Get mask for cells that are either null or nan
    bad_mask = np.logical_or(pd.isnull(df),pd.isna(df))
    idx = list(df.index)
    final_mask = []
    for i in range(len(bad_mask)):

        # Grab row
        row = df.loc[idx[i],:]

        # Get all non-null and non-na entries from the row
        mask = np.logical_not(np.logical_or(pd.isnull(row),pd.isna(row)))
        not_null = row[mask]

        # Default to not keeping a row
        keep = False

        # If there are any non-null entries left...
        if len(not_null) > 0:

            # Loop will end up with keep = False if only columns remaining are
            # type string and "".
            for n in not_null:

                # If an empty string, this is basically a null.
                if type(n) is str and n.strip() == "":
                    continue

                # If we get here, we have some kind of non-empty cell. Break out
                # and keep row
                keep = True
                break

        final_mask.append(keep)

    final_mask = np.array(final_mask,dtype=bool)

    # Let user know we're dropping lines
    if np.sum(np.logical_not(final_mask)) > 0:
        lines = ",".join(["{}".format(i) for i in df.index[np.logical_not(final_mask)]])
        w = f"\nDropping apparently empty lines ({lines})\n\n"
        print(w)

    # Drop rows that are all empty
    df = df.loc[df.index[final_mask],:]

    # -------------------------------------------------------------------------
    # Process keep column

    # Try to extract keep column
    try:
        keep = df.loc[:,"keep"]
    # No column: create one and set all to True
    except KeyError:
        keep = None
        df["keep"] = np.ones(len(df),dtype=bool)

    if keep is not None:

        # Do a pass trying to infer the datatype of keep. (This is useful if we
        # dropped empty rows that made the original pandas read this column in
        # as a mix of bool and object).
        df.loc[:,"keep"] = df.keep.infer_objects()

        # If it's not a boolean column, try to turn into one
        if not np.dtype(keep.dtypes) is np.dtype(bool):

            # Base message. If everything works great, let user know what
            # happened as warning. If things go awry, use as start of error
            # message
            w = "\n\n"
            w += "The 'keep' column must be boolean (True/False). pandas did\n"
            w += "not recognize the column as boolean, so we're parsing it\n"
            w += "manually by looking for 0/1, yes/no, true/false, etc.\n\n"

            new_keep = []
            look_for_true = re.compile("[1yt]",re.IGNORECASE)
            look_for_false = re.compile("[0nf]",re.IGNORECASE)
            for k in df.keep:
                if type(k) is bool:
                    is_true = True and k
                    is_false = not is_true
                    looks_like_a = "bool"
                elif type(k) is str:
                    is_true = look_for_true.search(k) is not None
                    is_false = look_for_false.search(k) is not None
                    looks_like_a = "string"
                elif type(k) is int:
                    is_true = (k != 0)
                    is_false = (k == 0)
                    looks_like_a = "int"
                elif type(k) is float:
                    is_true = np.logical_not(np.isclose(k,0))
                    is_false = np.isclose(k,0)
                    looks_like_a = "float"
                else:
                    w += f"Could not figure out how to parse '{k}'\n\n"
                    raise ValueError(w)

                if (is_true and is_false) or (not is_true and not is_false):
                    w += f"Trying to parse '{k}' as a {looks_like_a}, but\n"
                    w += f"could not figure out whether true or false.\n\n"
                    raise ValueError(w)
                else:
                    new_keep.append(is_true)

            # Record newly boolean-ized values
            df.loc[:,"keep"] = np.array(new_keep,dtype=bool)

            # Let user know we manually parsed the keep column...
            print(w)

    # -------------------------------------------------------------------------
    # Process uid column

    warn_uid = ""

    # Make sure the dataframe has a uid column
    try:
        uid = df.loc[:,"uid"]
    except KeyError:
        uid =  None
        warn_uid += "  + No 'uid' column found. Creating column.\n"
        df["uid"] = _private.generate_uid(len(df))

    if uid is not None:

        # Look for a-z and A-Z (only allowed characters in uid)
        p = re.compile("[^a-z]",re.IGNORECASE)

        final_uid = []

        # Go through uid...
        for u in uid:

            # Force to be a spring
            if type(u) is not str:
                u = str(u)

            # Disallowed character(s) in uid or uid not right length
            if p.search(u) or len(u) != 10:

                new_uid = _private.generate_uid()
                warn_uid += f"  + uid '{u}' is invalid. Replacing with '{new_uid}'\n"
                final_uid.append(new_uid)

            else:
                final_uid.append(u)

        df.loc[:,"uid"] = final_uid

    # Make sure uid are unique
    uid, counts = np.unique(df.loc[:,"uid"],return_counts=True)
    for u in uid[counts != 1]:
        warn_uid = f"  + uid '{u}' is not unique. This will be replaced with a\n"
        warn_uid += "    new unique uid for each sequence.\n"

        mask = df.loc[:,"uid"] == u
        df.loc[mask,"uid"] = _private.generate_uid(sum(mask))

    if warn_uid != "":
        w = "\n\n"
        w += "'uid' column was invalid. topiary has fixed the problems noted below.\n"
        w += "If you have already generated phylogenetic trees using a previous\n"
        w += "version of this dataframe, **those trees are no longer compatible\n"
        w += "with this dataframe**. To ensure compatibility, fix the problems\n"
        w += "noted below and then re-read the dataframe into topiary. If you\n"
        w += "have not already generated trees (or plan to generate new ones from\n"
        w += "scratch) you may safely disregard this warning.\n"
        w += "\n"
        w += warn_uid
        print(w)

    # -------------------------------------------------------------------------
    # Check format of ott column if present

    try:
        ott = df.loc[:,"ott"]
    except KeyError:
        ott = None

    if ott is not None:

        for index in df.index:

            # Get ott id
            o = df.loc[index,"ott"]

            # Missing value is okay --> make sure it's a pandas NULL datatype
            if pd.isna(o) or pd.isnull(o) or o == "":
                df.loc[index,"ott"] = pd.NA
                continue

            # If there is an ott that is not a null, make sure it is sane and
            # reasable.
            failed = False
            if type(o) is str and o[:3] == "ott":
                try:
                    int(o[3:])
                except ValueError:
                    failed = True
            else:
                failed = True

            if failed:
                err = "\n\nott column must have format 'ottINTEGER', where \n"
                err += "INTEGER is a the integer OTT accession number.\n"
                raise ValueError(err)


    # Make sure columns have order nickname, keep, species, name, sequence
    columns = list(df.columns)
    if "nickname" in columns:
        output_order = ["nickname"]
    else:
        output_order = []
    output_order.extend(["keep","species","name","sequence"])

    for c in columns:
        if c not in output_order:
            output_order.append(c)

    return df.loc[:,output_order]

def create_nicknames(df,
                     aliases=None,
                     source_column="name",
                     output_column="nickname",
                     overwrite_output=False,
                     auto_remove_common=True):
    """
    Create a nickname column that has a clean nickname for each sequence,
    generated by parsing "name" column.

    Parameters
    ----------
        df: topiary dataframe
        aliases: dictionary for creating standardized nicknames from input names.
                 Key specifies what should be output, values degenerate names
                 that map back to that key.  For example:
                     "S100A9":("S100-A9","S100 A9","S-100 A9")
                 would replace "S100-A9", "S100 A9", and "S-100 A9" with "S100A9"
        source_column: source_column in dataframe to use to generate a nickname
        output_column: place newly constructed nickname in this column.
        overwrite_output: overwrite an existing output column
        auto_remove_common: automatically remove some common words and replace
                            with "". Words to remove are: protein, product,
                            predicted, 'crystal structure of', hypothetical.

    Returns
    -------
        Copy of dataframe with new nickname column
    """

    # Working on a copy of the dataframe
    df = df.copy()

    # Check to make sure the specified source column exists
    try:
        df.loc[:,source_column]
    except KeyError:
        err = f"\ndataframe does not have source_column '{source_column}'\n\n"
        raise ValueError(err)

    # Make sure the output_column does not exist
    if output_column in _private.reserved_columns:
        err = f"\n'{output_column}' is a reserved column name. Please choose\n"
        err += "another column name.\n\n"
        raise ValueError(err)

    try:
        df.loc[:,output_column]
        if not overwrite_output:
            err = f"\ndataframe already has output_column '{output_column}'.\n"
            err += "To overwrite set overwrite_output = True\n\n"
            raise ValueError(err)
    except KeyError:
        pass

    # If aliases were not defined, set to an empty dictionary
    if aliases is None:
        aliases = {}

    # Check aliases data type
    if type(aliases) is not dict:
        err = "\naliases must be a dictionary\n\n"
        raise ValueError(err)

    # Work on a copy of the aliases dictionary
    aliases = copy.deepcopy(aliases)

    # Go through all alias keys
    for k in aliases:

        # Make sure the key is a string
        if type(k) is not str:
            err = "\naliases must be a dictionary keying strings to tuples\n"
            err += f"of strings. Alias '{k}' not recognized.\n\n"
            raise ValueError(err)

        # Make sure value is an iterable
        if hasattr(aliases[k],"__iter__"):

            # Wrap a naked string in a tuple, turning it into an iterable
            if type(aliases[k]) is str:
                aliases[k] = (aliases[k],)

            # Go through every alias and make sure it is a list of strings
            for a in aliases[k]:
                if type(a) is not str:
                    err = "\naliases must be a dictionary keying strings to tuples\n"
                    err += f"of strings. Alias '{a}' not recognized.\n\n"
                    raise ValueError(err)

            # Convert aliases values to tuples
            aliases[k] = tuple(aliases[k])


        else:
            err = "\naliases must be a dictionary keying strings to tuples\n"
            err += f"of strings. Value '{aliases[k]}' not recognized.\n\n"
            raise ValueError(err)


    # If automatically removing words, add those...
    if auto_remove_common:

        try:
            annihilate_words = list(aliases[""])
        except KeyError:
            annihilate_words = []

        annihilate_words.extend(["protein",
                                 "product",
                                 "PREDICTED:",
                                 "Crystal structure of",
                                 "hypothetical"])

        aliases[""] = tuple(annihilate_words)

    # Get entries from source column
    output = [entry for entry in df.loc[:,source_column]]

    # For every alias
    for a in aliases:

        # For every pattern corresponding to that alias
        for p in aliases[a]:

            # Replace patterns with alias for every line
            ap = re.compile(p,re.IGNORECASE)
            for i in range(len(output)):
                output[i] = ap.sub(a,output[i])

    # Clean up double spaces and leading/trailing spaces
    dbl = re.compile("  ")
    for i in range(len(output)):
        output[i] = dbl.sub(" ",output[i]).strip()

    df.loc[:,output_column] = output

    return check_topiary_dataframe(df)

def get_ott_id(df,phylo_context="All life"):
    """
    Get open taxonomy of life ids for all species in data frame.

    phylo_context: string. used to limit species seach for looking up species
                   ids on open tree of life.  To get latest strings recognized
                   by the database, use the following code:

                   ```
                   from opentree import OT
                   print(OT.tnrs_contexts().response_dict)
                   ```

                   As of 2021-08-16, the following are recognized. You can use
                   either the keys or values in this dictionary.

                   {'ANIMALS': ['Animals','Birds','Tetrapods','Mammals',
                                'Amphibians','Vertebrates','Arthropods',
                                'Molluscs','Nematodes','Platyhelminthes',
                                'Annelids','Cnidarians','Arachnids','Insects'],
                    'FUNGI': ['Fungi', 'Basidiomycetes', 'Ascomycetes'],
                    'LIFE': ['All life'],
                    'MICROBES': ['Bacteria','SAR group','Archaea','Excavata',
                                 'Amoebozoa','Centrohelida','Haptophyta',
                                 'Apusozoa','Diatoms','Ciliates','Forams'],
                    'PLANTS': ['Land plants','Hornworts','Mosses','Liverworts',
                               'Vascular plants','Club mosses','Ferns',
                               'Seed plants','Flowering plants','Monocots',
                               'Eudicots','Rosids','Asterids','Asterales',
                               'Asteraceae','Aster','Symphyotrichum',
                               'Campanulaceae','Lobelia']}
    """

    # Make sure this is a topiary dataframe
    df = check_topiary_dataframe(df)

    # Get ott id using open try of life
    df = opentree.get_ott_id(df,context_name=phylo_context)

    return df
