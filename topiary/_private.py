
import pandas as pd
import numpy as np

import re, sys, os, string, random, pickle, io, urllib, http

required_columns = ["species","name","sequence"]
reserved_columns = required_columns[:]
reserved_columns.extend(["uid","ott","alignment"])

def generate_uid(number=1):
    """
    Generate a unique uid. This will be a 10 character random combination of
    ascii letters.

    number: number of uid to generate. if 1, return a single uid. if > 1,
            return a list of uid.
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

def to_pretty(row):
    """
    Given a pandas Series, create pretty output.
    """

    try:
        pretty = f"{row.uid}|{row.nickname}|{row.species}"
    except AttributeError:
        try:
            pretty = f"{row.uid}|{row.species}"
        except AttributeError:
            err = "\n\nrow does not have all required attributes:"
            err += " (uid, species)\n"
            raise ValueError(err)

    pretty = re.sub("[,:;\"\'\(\)\.]","-",pretty)

    return pretty
