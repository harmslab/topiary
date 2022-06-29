"""
Construct a command line argument parser for a function, parse the command
line arguments, and then run the function.
"""

import sys, inspect, argparse, re, os

class IterFromFile(argparse.Action):
    """
    Argparse action for intelligently handling list-like arguments. Processes
    arguments assuming either 1) multiple arguments that are then put into a
    list of appropriate type or 2) a single file that has values on each line.
    When reading the file, class will delete everything after "#" on a line and
    will skip blank lines.
    """
    def __init__(self, option_strings, dest, value_type, **kwargs):
        # The "value_type" argument is the value that every entry in the iterable
        # will be cast to.
        self._value_type = value_type
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):

        # Name of type (i.e. int, str, etc.)
        type_name = self._value_type.__name__

        # If there is a single argument that is a file, get the values from that
        # file.
        read_from_file = False
        if len(values) == 1:
            if os.path.isfile(values[0]):
                read_from_file = True

        # If reading from a file...
        if read_from_file:

            final_values = []
            with open(values[0]) as f:
                for line in f:

                    # Split on # for comments
                    v = line.split("#")[0].strip()

                    # If line is blank, skip
                    if v == "": continue

                    try:
                        final_values.append(self._value_type(v))
                    except (TypeError,ValueError):

                        err = f"\n\Reading '{option_string} {values[0]}' as a file.\n"
                        err += f"Could not parse line '{v}' as '{type_name}'. {option_string} should\n"
                        err += f"either have a single argument that points to a file with one '{type_name}'\n"
                        err += f"on each line or have one or more arguments that can be interpreted as\n"
                        err += f"{type_name}. Examples:\n\n"
                        err += f"    {option_string} SOME_FILE\n"
                        err += f"    {option_string} {type_name}1 {type_name}2 ...\n\n"
                        raise ValueError(err)

        # If parsing arguments directly from the command line
        else:

            final_values = []
            for v in values:
                try:
                    final_values.append(self._value_type(v))
                except (TypeError,ValueError):

                    err = f"\n\nCould not parse argument '{v}' as '{type_name}'. {option_string} should\n"
                    err += f"either have a single argument that points to a file with one '{type_name}'\n"
                    err += f"on each line or have one or more arguments that can be interpreted as\n"
                    err += f"{type_name}. Examples:\n\n"
                    err += f"    {option_string} SOME_FILE\n"
                    err += f"    {option_string} {type_name}1 {type_name}2 ...\n\n"
                    raise ValueError(err)

        # Set the destination attribute in the name space to our parsed values.
        setattr(namespace, self.dest, final_values)

def wrap_function(fcn,argv=None,optional_arg_types={}):
    """
    Construct a command line argument parser for a function, parse the command
    line arguments, and then run the function.

    Parameters
    ----------
    fcn : function
        function to run.
    argv : list, optional
        arguments to parse. if None, use sys.argv[1:]
    optional_arg_types : dict, optional
        dictionary of arg types for arguments with None as their default in the
        function. If an argument where default is None is not in
        optional_arg_types, treat argument as str.

    Return
    ------
    argparse.ArgumentParser().parse_args(argv) : argparse.Namespace
        namespace shows how command line was parsed
    """

    # Get command line arguments
    if argv is None:
        argv = sys.argv[1:]

    # Get program name
    prog = re.sub("_","-",fcn.__name__)
    prog = f"topiary-{prog}"

    # Get description
    description = dict(inspect.getmembers(fcn))["__doc__"]
    description = re.sub(":code:","",description)

    # Build parser
    parser = argparse.ArgumentParser(prog=prog,
                                     description=description,
                                     formatter_class=argparse.RawTextHelpFormatter)

    # Build parser arguments using signature of fcn
    param = inspect.signature(fcn).parameters
    for p in param:

        # If no default specified, make required and move on.
        if param[p].default is param[p].empty:
            parser.add_argument(p)
            continue

        # For default is None args, parse as optional_arg_types or str
        if param[p].default is None:
            try:
                arg_type = optional_arg_types[p]
            except KeyError:
                arg_type = str

        # Otherwise, just grab the type
        else:
            arg_type = type(param[p].default)

        # bool
        kwargs = {}
        if arg_type is bool:
            if param[p].default is True:
                kwargs["action"] = "store_false"
            else:
                kwargs["action"] = "store_true"

        # non-str iterable. either read from a file or take + args in a row
        elif hasattr(arg_type,"__iter__") and arg_type is not str:
            kwargs["type"] = str
            kwargs["default"] = param[p].default
            kwargs["action"] = IterFromFile
            kwargs["value_type"] = type(param[p].default[0])
            kwargs["nargs"] = "+"

        # any other argument type
        else:
            kwargs["type"] = arg_type
            kwargs["default"] = param[p].default

        parser.add_argument(f"--{p}",**kwargs)

    # Parse stats
    args = parser.parse_args(argv)

    # Call function with kwargs
    try:
        fcn(**args.__dict__)
    except Exception as e:
        err = f"\n\nFunction {fcn.__name__} raised an error.\n\n"
        err += f"To see command line help, run {prog} --help\n\n"
        raise RuntimeError(err) from e

    return args
