__description__ = \
"""
Process paralog_patterns arguments and check their sanity.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-05-27"

import re, copy

def process_paralog_patterns(paralog_patterns,ignorecase=True,re_flags=None):
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
