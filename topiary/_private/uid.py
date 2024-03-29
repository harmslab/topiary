"""
Generate a unique uid. This will be a 10 character random combination of
ascii letters.
"""

import string, random

def generate_uid(number=1):
    """
    Generate a unique uid. This will be a 10 character random combination of
    ascii letters.

    Parameters
    ----------
    number : int, defulat=1
        number of uid to generate. if 1, return a single uid. if > 1, return a
        list of uid.

    Return
    ------
    out: str or list
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
