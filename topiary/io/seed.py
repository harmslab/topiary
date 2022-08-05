"""
Load a seed data frame and extract information needed to set up the
ancestral sequence reconstruction calculation.
"""

import topiary
from topiary._private import check
from topiary.external.opentree import species_to_ott, ott_resolvable

import numpy as np
import pandas as pd

import re, itertools
from string import ascii_lowercase, digits

compile_err = \
"""
Could not compile regular expression '{}'. Topiary uses regular expressions
constructed from user aliases to create human-readable sequence names and to
identify reciprocal blast hits. Topiary is struggling to create a valid, unique
regular expression for protein '{}'. Please make sure that none of the aliases
are the same between this protein and other proteins in the seed dataframe.
Also check that the aliases are not complete subsets of one another. For
example: one protein has the alias 'S100' and another 'S100A5'. Such alias
subsets usually work; however, because regular expression construction is
failing for this case, try avoiding such aliases (at least for troubleshooting).
Finally, try removing parentheses, backslashes and other special characters
from your aliases.
"""

def _get_alias_regex(some_string,spacers=[" ","-","_","."]):
    """
    Get regular expressions describing a string with variable spacers. For
    example LY96, LY-96, LY 96, LY_96. This function looks for elements in
    spacers or transitions between letters and numbers and builds regular
    expressions to match those patterns. It makes all letters lowercase
    and thus expects to be used in a regular expression ignoring case.

    Examples (for spacers = [" ","-","_","."]):
        "LY96" -> "ly[\ \-_\.]*96"
        "MD-2" -> "md[\ \-_\.]*2"
        "myeloid factor 2" -> "myeloid[\ \-_\.]*factor[\ \-_\.]*2"
        "CDR2L" -> "cdr[\ \-_\.]*2[\ \-_\.]*L"

    Parameters
    ----------
    some_string : str
        string to process
    spacers: list, default=[" ","-","_","."]
        list of characters to recognize as spacers

    Return
    ------
        regular expression as a string
    """

    # Regular expression describing spacers
    spacer_re = "".join([re.escape(s) for s in spacers])
    spacer_re = f"[{spacer_re}]*"

    # Make lowercase and strip leading/trailing white space
    lower = some_string.lower().strip()

    # Get digit, letter, and spacer characters
    is_digit = np.zeros(len(lower),dtype=bool)
    is_letter = np.zeros(len(lower),dtype=bool)
    is_spacer = np.zeros(len(lower),dtype=bool)
    for i, s in enumerate(lower):
        if s in digits:
            is_digit[i] = True
        elif s in ascii_lowercase:
            is_letter[i] = True
        elif s in spacers:
            is_spacer[i] = True
        else:
            pass

    # Go through sequence and build pattern character by character. Each block
    # separated by space will end up in its own list within the lists
    pattern = [[]]
    for i in range(len(lower)-1):

        # If we hit a spacer, put in non-greedy match
        if is_spacer[i]:
            pattern.append([])
            continue

        # Append the letter/digit we saw
        pattern[-1].append(lower[i])

        # Digit followed by letter
        if is_digit[i] and is_letter[i+1]:
            pattern.append([])

        # Letter followed by digit
        if is_letter[i] and is_digit[i+1]:
            pattern.append([])

    # Last letter
    pattern[-1].append(lower[-1])

    # Join letters within blocks and escape any regex symbols lurking in there
    pattern = [re.escape("".join(p)) for p in pattern if len(p) != 0]

    # Join pattern on spacer regular expression
    pattern = spacer_re.join(pattern)

    return pattern


def _build_alias_regex(alias_dict,spacers=[" ","-","_","."]):
    """
    Build regex to look for aliases when assigning proteins names from raw
    NCBI description strings and/or doing reciprocal blast.

    Parameters
    ----------
    alias_dict : dict
        dictionary keying protein names to aliases extracted from seed
        dataframe.
    spacers: list, default=[" ","-","_","."]
        list of characters to recognize as spacers

    Returns
    -------
    paralog_patterns : dict
        dictionary of compiled regular expressions to use to try to match
        paralogs. Keys are paralog names; values are regular expressions.
    """

    # List of unique, sorted names
    names = list(set(alias_dict.keys()))
    names.sort()

    # --------------------------------------------------------------------------
    # Create sorted list of unique alias names, including name itself as one of
    # the aliases.

    alias_regex = {}
    full_regex = {}
    for name in alias_dict:

        alias_dict[name].append(name.lower())
        alias_dict[name] = list(set(alias_dict[name]))
        alias_dict[name].sort()

        alias_regex[name] = []
        for alias in alias_dict[name]:
            alias_regex[name].append(_get_alias_regex(alias,spacers))

        # Get rid of any regex that are duplicated after this (ly96 vs ly-96,
        # for example)
        alias_regex[name] = list(set(alias_regex[name]))
        alias_regex[name].sort()

        full_regex[name] = "|".join(alias_regex[name])

    # --------------------------------------------------------------------------
    # Look for identical aliases for more than one proteins

    bad_overlap = []
    for i in range(len(alias_dict)):
        ni = names[i]
        ai = set(alias_dict[ni])
        for j in range(i+1,len(alias_dict)):
            nj = names[j]
            aj = set(alias_dict[nj])

            overlap = ai.intersection(aj)
            if len(overlap) > 0:
                for ovp in overlap:
                    bad_overlap.append((ni,nj,ovp))

    if len(bad_overlap) > 0:
        err = f"\n\nDifferent proteins have have identical aliases. (Note that\n"
        err += "aliases are case-insensitive). The following proteins have\n"
        err += "identical aliases:\n\n"
        for bad in bad_overlap:
            err += f"    '{bad[0]}' & '{bad[1]}' --> '{bad[2]}'\n"
        err += "\n\n"
        raise ValueError(err)

    # --------------------------------------------------------------------------
    # Add negative matching to prevent cross-matching. Check to see if the regex
    # for each protein picks up aliases from other protein(s). If so, put a
    # *negative* match  against that alias. For example, if the regex looks for
    # "ly96" for protein A it would match alias "ly96-2" for protein B. To
    # prevent this cross-match, this code will build a regex like this:
    # ^(?!.*(ly96-2)).*(ly96), which says "look for strings that do not match
    # ly96-2 but DO match ly96.

    neg_match_added = []
    negative_match = {}
    for name in names:

        negative_match[name] = []
        try:
            pattern = re.compile(full_regex[name],flags=re.IGNORECASE)
        except re.error:
            err = compile_err.format(full_regex[name],name)
            raise RuntimeError(err)

        for other_name in names:

            if name == other_name:
                continue

            for i in range(len(alias_regex[other_name])):

                other_alias = alias_dict[other_name][i]
                if pattern.search(other_alias):
                    negative_match[name].append(alias_regex[other_name][i])
                    neg_match_added.append((name,other_alias,other_name))

        negative_match[name] = list(set(negative_match[name]))
        negative_match[name].sort()
        negative_match[name] = "|".join(negative_match[name])

        if negative_match[name] != "":
            neg = negative_match[name]
            pos = full_regex[name]
            full_regex[name] = f"^(?!.*({neg})).*({pos})"

    if len(neg_match_added) > 0:

        w = ["\nRegular expression built from aliases from one protein match aliases"]
        w.append("from another protein. These are listed below. Topiary will now")
        w.append("modify the regular expression to attempt to prevent this cross-matching.\n")
        w.append("Cross-matches are:")
        for n in neg_match_added:
            w.append(f"The protein '{n[0]}' aliases match alias '{n[1]}' for protein '{n[2]}'")
        w.append("\n")
        w = "\n".join(w)

        print(w)

    # --------------------------------------------------------------------------
    # compile final regular expressions

    for name in names:
        try:
            full_regex[name] = re.compile(full_regex[name],flags=re.IGNORECASE)
        except re.error:
            err = compile_err.format(full_regex[name],name)
            raise RuntimeError(err)


    # --------------------------------------------------------------------------
    # Validate final regexes, making sure they match aliases for self but non-
    # self.

    correct_match = []
    incorrect_match = []
    correct_missing = []
    incorrect_missing = []

    for name in names:

        pattern = full_regex[name]
        for other_name in names:

            for other_alias in alias_dict[other_name]:

                matches = not (pattern.search(other_alias) is None)

                if matches:
                    if name == other_name:
                        correct_match.append((name,other_alias,other_name))
                    else:
                        incorrect_match.append((name,other_alias,other_name))

                else:
                    if name == other_name:
                        incorrect_missing.append((name,other_alias,other_name))
                    else:
                        correct_missing.append((name,other_alias,other_name))

    if len(incorrect_match) > 0 or len(incorrect_missing) > 0:
        err = "\nAliases are ambiguous and cannot be used to ..."

        if len(incorrect_match) > 0:
            err += "Regular expressions that match more than one protein:\n\n"
            for m in incorrect_match:
                err += f"    '{m[0]}' regex matches alias '{m[1]}' for protein '{m[2]}'\n"
            err += "\n\n"

        if len(incorrect_missing) > 0:
            err += "Regular expressions that do not match aliases for their own protein:\n\n"
            for m in incorrect_missing:
                err += f"    '{m[0]}' regex does not match alias '{m[1]}'\n"
            err += "\n\n"

        raise ValueError(err)

    if len(neg_match_added) > 0:
        print("Success. Topiary was able to modify the regular expressions to")
        print("prevent cross-matching between protein aliases.")
        print("",flush=True)

    return full_regex


def read_seed(df):
    """
    Read a seed data frame and extract alias patterns and key species.

    Parameters
    ----------
    df: pandas.DataFrame or str
        seed dataframe containing seed sequences to launch the analysis. df can
        be a pandas dataframe or a string pointing to a spreadsheet file.

    Returns
    -------
    topiary_dataframe : pandas.DataFrame
        new topiary dataframe built from the seed dataframe
    key_species : numpy.array
        list if key species to keep during the analysis
    paralog_patterns : list
        list of compiled regular expressions to use to try to match paralogs.

    Notes
    -----

    The seed dataframe is expected to have at least four columns:

    + :code:`species`: species names for seed sequences in binomial format (i.e.
      Homo sapiens or Mus musculus)
    + :code:`name`: name of each sequence (i.e. LY96)
    + :code:`aliases`: other names for this sequence found in different
      databases/species, separated by ; (i.e. LY96;MD2;ESOP1)
    + :code:`sequence`: amino acid sequences for these proteins.

    It may have one other optional column:

    + :code:`recip_blast`: True/False. Indicates whether or not this species
      should be used as a key species for reciprocal BLASTing.

    Other columns in the dataframe are kept but not used by topiary.
    """

    # -----------------------------------------------------------------------
    # Load dataframe and check it's sanity

    if type(df) is not type(pd.DataFrame()):

        if type(df) is str:
            extension = df.split(".")[-1].strip().lower()
            if extension in ["xls","xlsx"]:
                df = pd.read_excel(df)
            elif extension == "csv":
                df = pd.read_csv(df,sep=",")
            elif extension == "tsv":
                df = pd.read_csv(df,sep="\t")
            else:
                df = pd.read_csv(df,sep=None,engine="python")
        else:
            err = f"Could not figure out how to read df '{df}'\n\n"
            raise ValueError(err)

    required_columns = ["species","name","aliases","sequence"]
    df.columns = [c.lower().strip() for c in df.columns]
    for c in required_columns:

        try:
            df[c]
        except KeyError:
            err = f"\nDataframe must have a {c} column. Required columns are:\n"
            for r in required_columns:
                err += f"    {r}\n"
            err += "\n"
            raise ValueError(err)

        # Strip extra leading/trailing white space
        df.loc[:,c] = df.loc[:,c].str.strip()


    # Make sure all input species are found in the OTT database
    bad_species = []
    ott_list, species_list, _ = species_to_ott(np.unique(df.loc[:,"species"]))
    for i in range(len(species_list)):
        if ott_list[i] is None:
            bad_species.append(species_list[i])

    if len(bad_species) > 0:
        err = "\nNot all input species were found in the Open Tree of Life\n"
        err += "database. To troubleshoot the problem, you can visit\n"
        err += "https://tree.opentreeoflife.org/taxonomy/browse and search for\n"
        err += "the following species manually:\n\n"
        for b in bad_species:
            err += f"    '{b}'\n"
        raise ValueError(err)

    # Make sure all input species can be resolved on the OTT synthetic tree
    resolved = ott_resolvable(ott_list)
    unresolvable_species = []
    for i in range(len(resolved)):
        if not resolved[i]:
            unresolvable_species.append((species_list[i],ott_list[i]))

    if len(unresolvable_species) > 0:
        err = "\nNot all input species can be resolved on the latest Open Tree of\n"
        err += "Life synthetic tree. To troubleshoot the problem, you can visit\n"
        err += "https://tree.opentreeoflife.org/taxonomy/browse and search for\n"
        err += "the OTT accession numbers of these species:\n"

        for b in unresolvable_species:
            err += f"    '{b[0]}' (OTT '{b[1]}')\n"
        raise ValueError(err)

    # -----------------------------------------------------------------------
    # Get key_species

    # Look for key_species. If present, make sure it is bool. If not, add it
    # as all True
    if "key_species" in df.columns:
        df.loc[:,"key_species"] = check.column_to_bool(df.loc[:,"key_species"],
                                                       "key_species")
    else:
        df["key_species"] = True

    key_species = np.unique(df.loc[df.loc[:,"key_species"],"species"])
    key_species.sort()

    # -----------------------------------------------------------------------
    # Get regular expression for searching for paralogs

    alias_dict = {}
    for idx in df.index:

        # Construct list of all aliases given by the user for a given name
        key = str(df.loc[idx,"name"]).strip()
        values = [v.strip().lower() for v in df.loc[idx,"aliases"].split(";")]

        # Get rid of ""
        values = [v for v in values if v != ""]

        try:
            alias_dict[key].extend(values)
        except KeyError:
            alias_dict[key] = values[:]

    paralog_patterns = _build_alias_regex(alias_dict)

    # Load newly created paralog patterns into aliases column
    new_aliases = []
    for idx in df.index:

        name = str(df.loc[idx,"name"]).strip()
        new_aliases.append(paralog_patterns[name].pattern)

    df.loc[:,"aliases_regex"] = new_aliases


    # -----------------------------------------------------------------------
    # Convert dataframe into a topiary dataframe

    # add always keep column to this dataframe
    df["always_keep"] = True

    # Drop on uid column to avoid warning when topiary creates in arg_processor
    try:
        df["uid"]
    except KeyError:
        df["uid"] = topiary._private.generate_uid(len(df["always_keep"]))

    df = check.check_topiary_dataframe(df)

    return df, key_species, paralog_patterns
