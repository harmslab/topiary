"""
Convert a paralog patterns dictionary with lists of patterns as values into a
dictionary with regex as values. 
"""

import numpy as np

import re
import copy
from string import ascii_letters, digits

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

def _get_alias_regex(some_string,
                     spacers=[" ","-","_","."],
                     ignorecase=True):
    """
    Get regular expressions describing a string with variable spacers. For
    example LY96, LY-96, LY 96, LY_96. This function looks for elements in
    spacers or transitions between letters and numbers and builds regular
    expressions to match those patterns. 

    Examples (for spacers = [" ","-","_","."]):
        "LY96" -> "ly[\\ \\-_\\.]*96"
        "MD-2" -> "md[\\ \\-_\\.]*2"
        "myeloid factor 2" -> "myeloid[\\ \\-_\\.]*factor[\\ \\-_\\.]*2"
        "CDR2L" -> "cdr[\\ \\-_\\.]*2[\\ \\-_\\.]*L"

    Parameters
    ----------
    some_string : str
        string to process
    spacers: list, default=[" ","-","_","."]
        list of characters to recognize as spacers
    ignorecase : bool, default=True
        when compiling regex, whether or not to ignore case

    Return
    ------
        regular expression as a string
    """

    # Regular expression describing spacers
    spacer_re = "".join([re.escape(s) for s in spacers])
    spacer_re = f"[{spacer_re}]*"

    # Make lowercase and strip leading/trailing white space
    if ignorecase:
        this_string = some_string.lower().strip()
    else:
        this_string = some_string.strip()

    # Get digit, letter, and spacer characters
    is_digit = np.zeros(len(this_string),dtype=bool)
    is_letter = np.zeros(len(this_string),dtype=bool)
    is_spacer = np.zeros(len(this_string),dtype=bool)
    for i, s in enumerate(this_string):
        if s in digits:
            is_digit[i] = True
        elif s in ascii_letters:
            is_letter[i] = True
        elif s in spacers:
            is_spacer[i] = True
        else:
            pass

    # Go through sequence and build pattern character by character. Each block
    # separated by space will end up in its own list within the lists
    pattern = [[]]
    for i in range(len(this_string)-1):

        # If we hit a spacer, put in non-greedy match
        if is_spacer[i]:
            pattern.append([])
            continue

        # Append the letter/digit we saw
        pattern[-1].append(this_string[i])

        # Digit followed by letter
        if is_digit[i] and is_letter[i+1]:
            pattern.append([])

        # Letter followed by digit
        if is_letter[i] and is_digit[i+1]:
            pattern.append([])

    # Last letter
    pattern[-1].append(this_string[-1])

    # Join letters within blocks and escape any regex symbols lurking in there
    pattern = [re.escape("".join(p)) for p in pattern if len(p) != 0]

    # Join pattern on spacer regular expression
    pattern = spacer_re.join(pattern)

    return pattern


def _build_alias_regex(alias_dict,
                       spacers=[" ","-","_","."],
                       re_flags=None,
                       ignorecase=True):
    """
    Build regex to look for aliases when assigning proteins names from raw
    NCBI description strings and/or doing reciprocal blast. The alias_dict
    can either have a list of strings as values or a pre-compiled regular 
    expression for each value. If a pre-compiled expression, that expression is
    used as-is. Otherwise, the list of strings is compiled into a regular 
    expression. 

    Parameters
    ----------
    alias_dict : dict
        dictionary keying protein names to either a list of aliases (as strings)
        OR a pre-compiled regular expression. 
    spacers: list, default=[" ","-","_","."]
        list of characters to recognize as spacers
    re_flags : list, optional
        regular expression flags to pass to compile. None or list of
        of flags. Note, "ignorecase" takes precedence over re_flags.
    ignorecase : bool, default=True
        when compiling regex, whether or not to ignore case

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

        # If regex passed in, just use as-is
        if issubclass(type(alias_dict[name]),re.Pattern):
            alias_regex[name] = [alias_dict[name].pattern]
            full_regex[name] = alias_dict[name]

        # Otherwise, construct new regex
        else:

            if ignorecase:
                alias_dict[name].append(name.lower())
            else:
                alias_dict[name].append(name)

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
    # Look for identical aliases for more than one proteins. (Note this will
    # not catch identical pre-compiled regex). 

    bad_overlap = []
    for i in range(len(alias_dict)):
        ni = names[i]
        if issubclass(type(alias_dict[ni]),re.Pattern):
            continue

        ai = set(alias_dict[ni])
        for j in range(i+1,len(alias_dict)):
            nj = names[j]
            if issubclass(type(alias_dict[nj]),re.Pattern):
                continue
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
        if issubclass(type(full_regex[name]),re.Pattern):
            pattern = full_regex[name]
        else:         
            try:
                pattern = re.compile(full_regex[name],**re_kwargs)
            except re.error:
                err = compile_err.format(full_regex[name],name)
                raise RuntimeError(err)

        for other_name in names:

            if name == other_name:
                continue

            for i in range(len(alias_regex[other_name])):

                if issubclass(type(alias_dict[other_name]),re.Pattern):
                    other_alias = alias_dict[other_name].pattern
                else:
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

        if not issubclass(type(full_regex[name]),re.Pattern):
            try:
                full_regex[name] = re.compile(full_regex[name],**re_kwargs)
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

        if issubclass(type(alias_dict[name]),re.Pattern):
            continue

        pattern = full_regex[name]
        for other_name in names:

            if issubclass(type(alias_dict[other_name]),re.Pattern):
                continue

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


def _preprocess_alias_input(alias_dict):
    """
    Do error checking on an alias dictionary and turn into the expected format
    for _build_alias_regex.

    Parameters
    ----------
    alias_dict : dict
        dictionary to check

    Returns
    ------
    alias_dict : dict
        dictionary keying patterns to lists of strings 
    """

    # Make generic, informative, error when dealing with paralog_patterns
    generic_pp_error = ["alias_dict must be a dictionary keying paralog",
                        "names to patterns that match that paralog. For example,",
                        "the protein LY96 is annotated as LY96, MD2, MD-2,",
                        "ESOP1, Lymphocyte antigen 96, or myeloid differentiation 2.",
                        "A alias_dict dictionary for this protein would be",
                        "alias_dict = {'LY96:['MD2','MD-2','ESOP1',",
                        "                      'lymphocyte antigen 96',",
                        "                      'myeloid differentiation 2']}\n",
                        "The patterns can either be strings (i.e. 'MD-2') or",
                        "re.Pattern instances (i.e. re.compile('MD().*?)2'))."]
    generic_pp_error = "\n".join(generic_pp_error)

    # Check paralog_patterns data type
    if not issubclass(type(alias_dict),dict):
        err = f"\nalias_dict is not a dictionary\n\n{generic_pp_error}\n\n"
        raise ValueError(err)

    # Work on a copy of the paralog_patterns dictionary
    alias_dict = copy.deepcopy(alias_dict)

    # Go through all paralog_patterns keys
    for k in alias_dict:

        # Make sure the key is a string
        if not issubclass(type(k),str):
            err = f"\nparalog key '{k}' not recognized.\n\n{generic_pp_error}\n\n"
            raise ValueError(err)

         # If value is a compiled expression pattern, leave alone
        if issubclass(type(alias_dict[k]),re.Pattern):
            continue

        # If value is a string, turn into a list so we can iterate over it 
        if issubclass(type(alias_dict[k]),str):
            alias_dict[k] = [alias_dict[k]]
            continue

        # Make sure value is an iterable of strings or patterns
        if hasattr(alias_dict[k],"__iter__"):

            # Make sure there is at least one pattern
            if len(alias_dict[k]) == 0:
                err = f"\nparalog key '{k}' has no patterns.\n\n{generic_pp_error}\n\n"
                raise ValueError(err)

            alias_dict[k] = list(alias_dict[k])

            # Make sure every value in list is a string
            for i in range(len(alias_dict[k])):
                if issubclass(type(alias_dict[k][i]),str):
                    continue
                elif issubclass(type(alias_dict[k][i]),re.Pattern):
                    alias_dict[k][i] = alias_dict[k][i].pattern
                else:
                    err = f"\nalias_dict entry {alias_dict[k][i]} for key {k} is not\n"
                    err += "a str or re.Pattern\n\n"
                    err += f"{generic_pp_error}\n\n"
                    raise ValueError(err)

        # value is not iterable -- bad news
        else:
            err = f"\nvalue for alias_dict '{k}' should be list-like or a\n"
            err += "compiled regular expression.\n\n"
            err += f"{generic_pp_error}\n\n"
            raise ValueError(err)

    return alias_dict


def load_paralog_patterns(alias_dict,
                          spacers=[" ","-","_","."],
                          ignorecase=True,
                          re_flags=None):
    """
    Build regex to look for aliases when assigning protein names from raw
    NCBI description strings and/or doing reciprocal blast. The alias_dict
    can either have a list of strings as values or a pre-compiled regular 
    expression for each value. If a pre-compiled expression, that expression is
    used as-is. Otherwise, the list of strings is compiled into a regular 
    expression. 

    Parameters
    ----------
    alias_dict : dict
        dictionary keying protein names to either a list of aliases (as strings)
        OR a pre-compiled regular expression. 
    spacers: list, default=[" ","-","_","."]
        list of characters to recognize as spacers
    re_flags : list, optional
        regular expression flags to pass to compile. None or list of
        of flags. Note, "ignorecase" takes precedence over re_flags.
    ignorecase : bool, default=True
        when compiling regex, whether or not to ignore case

    Returns
    -------
    paralog_patterns : dict
        dictionary of compiled regular expressions to use to try to match
        paralogs. Keys are paralog names; values are regular expressions.
    """
    
    alias_dict = _preprocess_alias_input(alias_dict)

    pp = _build_alias_regex(alias_dict,
                            spacers=spacers,
                            ignorecase=ignorecase,
                            re_flags=re_flags)

    return pp

