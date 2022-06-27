"""
Private utility functions that are not publicly exposed in the API.
"""

import string, random, sys, inspect, argparse, re

required_columns = ["species","name","sequence"]
reserved_columns = required_columns[:]
reserved_columns.extend(["uid","ott","alignment","keep","always_keep"])

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

def wrap_function(fcn,argv=None,optional_arg_types={}):
    """
    A generalized main function that constructs a command line argument parser
    and then parses arguments for any function.

    fcn: function to run.
    argv: arguments to parse. if none, use sys.argv[1:]
    optional_arg_types: dictionary of arg types for arguments with None as
                        their default in the function. If an argument is
                        not in optional_arg_types, the default is to treat
                        argument as an int.
    """

    # Get command line arguments
    if argv is None:
        argv = sys.argv[1:]

    # Get program name
    prog = re.sub("_","-",fcn.__name__)
    prog = f"topiary-{prog}"

    # Get description
    description = dict(inspect.getmembers(fcn))["__doc__"]

    # Build parser
    parser = argparse.ArgumentParser(prog=prog,
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # Build parser arguments using signature of fcn
    param = inspect.signature(fcn).parameters
    for p in param:

        # If no default specified, make required
        if param[p].default is param[p].empty:
            parser.add_argument(p)

        # If default specified, make optional
        else:

            # For type == None args, parse as string
            if param[p].default is None:
                try:
                    arg_type = optional_arg_types[p]
                except KeyError:
                    arg_type = str

            else:
                arg_type = type(param[p].default)

            parser.add_argument(f"--{p}",
                                type=arg_type,
                                default=param[p].default)

    # Parse stats
    args = parser.parse_args(argv)

    # Call function with kwargs
    try:
        fcn(**args.__dict__)
    except Exception as e:
        err = f"\n\n{description}\n\nIncorrect arguments specified.\n\n"
        err += "See command line help and error raised above.\n\n"
        raise RuntimeError(err) from e
