"""
Function to check/process topiary dataframes.
"""

from topiary import _private
from .standard import column_to_bool

import pandas as pd
import numpy as np

import re

def check_topiary_dataframe(df):
    """
    Check to make sure topiary dataframe is valid, adding standard columns if
    necessary. Edits and returns a copy of the input dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe to validate

    Returns
    -------
    pandas.DataFrame
        validated copy of input dataframe

    Notes
    -----
    + Checks for required columns 'species', 'name', and 'sequence'. Casts these
      entries as strings and does not allow them to be None/NaN/pandas.NA.
    + Deletes duplicated rows.
    + Deletes empty rows.

    + `keep`:
        + Must be boolean.
        + If present but not boolean, try to coerce it into boolean.
            + if str, map regex [1yt] --> True and regex[0nf] --> False
            + if int, map !0 --> True and 0 --> False
            + if float, map !~0 --> True and ~0 --> False
        + If keep is not present, add and set all to True.
    + `uid`:
        + Must be unique str of length 10 consisting only of [a-z][A-Z]
        + If present:
            + If missing values, fill them in
            + If values do not meet rules, replace them with uid that do
            + If values are not unique, update them to be unique
        + If not present, create them.
    + `ott`: If present, makes sure it has the format ottINTEGER.
    + `alignment`: if present, makes sure that all columns are the same width,
      and removes any gaps-only columns.

    + Enforces following column order:
        + `nickname` (if present)
        + `keep`
        + `species`
        + `name`
        + `sequence`
        + all other columns, in order they came in
    """

    # Make sure type is right
    if type(df) is not pd.DataFrame:
        err = "\ndf must be a pandas dataframe.\n\n"
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
    # Process required columns. Enforce them as non-null and castable to string

    for c in _private.required_columns:

        try:
            df[c]
        except KeyError:
            err = f"\n\nA topiary dataframe must have a {c} column.\n\n"
            raise ValueError(err)

        # Do a pass trying to infer the datatype of c. (Important in case we
        # dropped empty rows)
        df.loc[:,c] = df.loc[:,c].infer_objects()
        casted = []
        for v in df.loc[:,c]:
            if pd.isna(v) or pd.isnull(v) or v is None or not v or v == "":
                err = f"\nMissing value in required column '{c}'.\n\n"
                raise ValueError(err)

            try:
                casted.append(str(v))
            except (TypeError,ValueError):
                err = f"\nColumn '{c}' must have only strings.\n\n"
                raise ValueError(err)

        df.loc[:,c] = casted

    # -------------------------------------------------------------------------
    # Process sequence column

    # stripping leading/trailing whitespace and removing any in the middle. This
    # is safe because operation above casted to string
    for idx in df.index:
        sequence = str(df.loc[idx,"sequence"]).strip()
        df.loc[idx,"sequence"] = "".join(sequence.split())

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

        # Validate this as a boolean column
        df.loc[:,"keep"] = column_to_bool(df.loc[:,"keep"],"keep")

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
