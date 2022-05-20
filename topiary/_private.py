
import string, random

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
