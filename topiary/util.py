
__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
utility functions for the topiary package
"""

from . import _private

import pandas as pd
import numpy as np

import re, sys, os, string, random, pickle, io, urllib, http, copy

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
    + 'alignment': if present, makes sure that all columns are the same width,
                   and removes any gaps-only columns.

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
        w = f"\nDropping empty lines ({lines})\n\n"
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
            if pd.isna(o) or pd.isnull(o) or o is None or not o or o in ["","none","None"]:
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

    # -------------------------------------------------------------------------
    # Check alignment column
    try:
        alignment = df.loc[:,"alignment"]
    except KeyError:
        alignment = None

    # If there is an alignment column...
    if alignment is not None:

        # Grab non-null alignment lines
        align_mask = np.logical_not(pd.isnull(df.alignment))
        align_mask = np.logical_and(align_mask,df.keep)

        # Make sure alignment column has at least one kept row
        if len(align_mask) > 0 and np.sum(align_mask) > 0:

            # Create matrix holding alignment, rows are sequences, columns are columns
            align_matrix = []
            for row in df.alignment[align_mask]:
                align_matrix.append(list(row))

            # Make sure all alignment rows are the same length
            unique_row_length = set([len(row) for row in align_matrix])
            if len(unique_row_length) != 1:
                err = "\nAll sequences in the 'alignment' column must have the\n"
                err += "same length\n\n"
                raise ValueError(err)

            # Convert to a matrix
            align_matrix = np.array(align_matrix)

            # Create mask for good columns -- columns with more than just "-"
            good_column_mask = []
            for i in range(align_matrix.shape[1]):
                u = np.unique(align_matrix[:,i])
                if len(u) == 1 and u[0] == "-":
                    good_column_mask.append(False)
                else:
                    good_column_mask.append(True)

            good_column_mask = np.array(good_column_mask,dtype=bool)

            # Whack out columns that are only "-"
            align_matrix = align_matrix[:,good_column_mask]

            # Convert alignment matrix back to an array of strings
            new_align = []
            for i in range(align_matrix.shape[0]):
                new_align.append("".join(align_matrix[i,:]))

            # Store in dataframe
            df.loc[align_mask,"alignment"] = new_align

    # -------------------------------------------------------------------------
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

def _compile_paralog_patterns(paralog_patterns,ignorecase=True,re_flags=None):
    """
    Internal function. Process a user-passed paralog patterns dictionary.

    Parameters
    ----------
        paralog_patterns: dictionary to check and compile
        ignorecase: when compiling regex, whether or not to ignore case
        re_flags: regular expression flags to pass to compile. None or list of
                  of flags. Note, "ignorecase" takes precedence over re_flags.

    Return
    ------
        list of patterns with format [(compiled_pattern,paralog),...]
    """

    # If nothing is passed in, return a list without throwing an error
    if paralog_patterns is None:
        return []

    # Make generic, informative, error when dealing with paralog_patterns
    generic_pp_error = ["paralog_patterns must be a dictionary keying paralog",
                        "names to patterns that match that paralog. For example,",
                        "the protein LY96 is annotated as LY96, MD2, MD-2,",
                        "ESOP1, Lymphocyte antigen 96, or myeloid differentiation 2.",
                        "A paralog_pattern dictionary for this protein would be",
                        "\nparalog_pattern = {'LY96:['MD2','MD-2','ESOP1',",
                        "                          'lymphocyte antigen 96',",
                        "                          'myeloid differentiation 2']}\n",
                        "The patterns can either be strings (i.e. 'MD-2') or",
                        "re.Pattern instances (i.e. re.compile('MD().*?)2')).",
                        "Note for advanced users, any strings passed in are",
                        "escaped. To use regex, pre-compile your regular",
                        "expressions."]
    generic_pp_error = "\n".join(generic_pp_error)

    # Deal with re_flags argument
    if re_flags is None:
        re_flags = []

    # Parse ignorecase argument
    if ignorecase:
        re_flags.append(re.IGNORECASE)

    # Assemble flags
    re_kwargs = {}
    if len(re_flags) > 0:
        re_kwargs["flags"] = re_flags[0]
        if len(re_flags) > 1:
            for f in re_flags[1:]:
                re_kwargs["flags"] = re_kwargs["flags"] | f

    # Check paralog_patterns data type
    if type(paralog_patterns) is not dict:
        err = "\nparalog_patterns is not a dictionary\n\n{generic_pp_error}\n\n"
        raise ValueError(err)

    # Work on a copy of the paralog_patterns dictionary
    pp = copy.deepcopy(paralog_patterns)

    # Go through all paralog_patterns keys
    patterns = []
    for k in paralog_patterns:

        # Make sure the key is a string
        if type(k) is not str:
            err = f"\nparalog key '{k}' not recognized.\n\n{generic_pp_error}\n\n"
            raise ValueError(err)

        # If value is a naked regular expression pattern or string, wrap it in a
        # list so the code can iterate over it.
        if type(pp[k]) in [re.Pattern,str]:
            pp[k] = [pp[k]]

        # Make sure value is an iterable
        if hasattr(pp[k],"__iter__"):

            # Make sure there is at least one pattern
            if len(pp[k]) == 0:
                err = f"\nparalog key '{k}' has no patterns.\n\n{generic_pp_error}\n\n"
                raise ValueError(err)

            # Compile patterns to search for in the source column
            to_compile = []
            for a in pp[k]:

                if type(a) is str:
                    to_compile.append(re.escape(a))
                elif type(a) is re.Pattern:
                    # NOTE: if user passes multiple compiled patterns with
                    # different flags, this will only use the flags for the
                    # last pattern on *all* regex passed in.
                    to_compile.append(a.pattern)
                    re_kwargs["flags"] = a.flags
                else:
                    err = f"\npattern '{a}' not recognized.\n\n{generic_pp_error}\n\n"
                    raise ValueError(err)

            patterns.append((re.compile("|".join(to_compile),**re_kwargs),k))

        # value is not iterable -- bad news
        else:
            err = f"\nvalue for paralog_pattern '{k}' should be list-like.\n\n"
            err += "{generic_pp_error}\n\n"
            raise ValueError(err)

    return patterns

def create_nicknames(df,
                     paralog_patterns,
                     source_column="name",
                     output_column="nickname",
                     separator="/",
                     unassigned_name="unassigned",
                     overwrite_output=False,
                     ignorecase=True):
    """
    Create a nickname column that has a friendly nickname for each sequence,
    generated by looking for patterns defined in the 'paralog_patterns'
    dictionary in 'source_column' column from the dataframe.

    Parameters
    ----------
        df: topiary dataframe
        paralog_patterns: dictionary for creating standardized nicknames from
                          input names. Key specifies what should be output, values
                          the a list of patterns that map back to that key.  For
                          example:

                          {"S100A9":["S100-A9","S100 A9","S-100 A9","MRP14"],
                           "S100A8":["S100-A8","S100 A8","S-100 A9","MRP8"]}

                          would assign "S100A9" to any sequence matching patterns
                          only from its list; "S100A8" to any sequence matching
                          patterns only from its list; and S100A9/S100A8 to any
                          sequence matching patterns from both lists.
        source_column: source column in dataframe to use to generate a nickname
                       (accessed via df.loc[:,source_column])
        output_column: column in which to store newly constructed nicknames
                       (accessed via df.loc[:,output_column])
        separator: "/" character to place between nicknames if more than one
                   pattern matches.
        unassigned_name: nickname to give sequences that do not match any of
                         the patterns.
        overwrite_output: boolean (default False). overwrite an existing output
                          column
        ignorecase: boolean (default True). Whether or not to ignore the case
                    of matches when assigning the nickname.

    Returns
    -------
        Copy of dataframe with new nickname column
    """

    # Parse 'df' argument
    if type(df) is not pd.DataFrame:
        err = "\ndf should be a pandas dataframe\n\n"
        raise ValueError(err)

    # Check to make sure the specified source column exists
    try:
        df.loc[:,source_column]
    except KeyError:
        err = f"\ndataframe does not have source_column '{source_column}'\n\n"
        raise ValueError(err)

    # Make sure the output column is not reserved by topiary
    if output_column in _private.reserved_columns:
        err = f"\n'{output_column}' is a reserved column name. Please choose\n"
        err += "another column name.\n\n"
        raise ValueError(err)

    # Make sure the output_column does not exist or that we're allowed to
    # overwrite
    try:
        df.loc[:,output_column]
        if not overwrite_output:
            err = f"\ndataframe already has output_column '{output_column}'.\n"
            err += "To overwrite set overwrite_output = True\n\n"
            raise ValueError(err)
    except KeyError:
        pass

    # check unassigned_name argument
    if type(unassigned_name) is not str:
        err = "\nunassigned_name should be a string.\n\n"
        raise ValueError(err)

    # check separator argument
    if type(separator) is not str:
        err = "\nseparator should be a string.\n\n"
        raise ValueError(err)

    patterns = _compile_paralog_patterns(paralog_patterns,ignorecase=ignorecase)

    # Get entries from source column
    source = [entry for entry in df.loc[:,source_column]]
    out = []
    for i, s in enumerate(source):

        out.append([])

        # Go over every patterns
        for p in patterns:

            # If we find the pattern in the entry, break; we only need
            # to note we found it once
            if p[0].search(source[i]):
                out[-1].append(p[1])
                break

        # Join list of hits
        out[-1] = f"{separator}".join(out[-1])

        # No hit, give unassigned_name
        if out[-1] == "":
            out[-1] = unassigned_name

    # Return an edited copy of the dataframe
    df = df.copy()
    df.loc[:,output_column] = out

    # Validate topiary dataframe to make sure not mangled; will also update
    # column order so nickname is early and thus in a user-friendly place

    return check_topiary_dataframe(df)
