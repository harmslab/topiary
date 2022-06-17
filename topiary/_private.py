__description__ = \
"""
Private utility functions that are not publicly exposed in the API.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-07"

import string, random

required_columns = ["species","name","sequence"]
reserved_columns = required_columns[:]
reserved_columns.extend(["uid","ott","alignment","length","keep","always_keep"])

# Data going into a newick tree can't have any of these symbols. We also reserve
# the '#' character for comments.
reserved_characters = ["(",")",";","#",":",",","'","\""]

def generate_uid(number=1):
    """
    Generate a unique uid. This will be a 10 character random combination of
    ascii letters.

    Parameters
    ----------
        number: number of uid to generate. if 1, return a single uid. if > 1,
                return a list of uid.

    Return
    ------
        uid or list of uid
    """

    if number < 1:
        err = "number must be 1 or more\n"
        raise ValueError(err)

    out = []
    for n in range(number):
        out.append("".join([random.choice(string.ascii_letters)
                            for _ in range(10)]))

    if number == 1:
        return out[0]

    return out

def create_pipeline_dict():
    """
    Create a dictionary with keys for column names and lists for values.
    This can be populated by a loop through sequences, followed by
    pandas.DataFrame(this_dictionary) to create a pandas data frame of
    the sort expected by the functions used in this module.

    Return
    ------
        dictionary with stereotyped keys
    """

    # Column names for output dictionary
    key_list = ["name","species","sequence","uid","keep"]

    return dict([(k,[]) for k in key_list])
