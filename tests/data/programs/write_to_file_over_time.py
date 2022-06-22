#!/usr/bin/env python
import sys, argparse, inspect, time

def time_waster(output_file,num_steps=5,time_step=0.2):
    """
    Function that that writes to output_file num_steps times, waiting time_step
    between each write.

    output_file: output file string
    num_steps: how many steps to take
    time_step: how long to wait between steps.
    """

    total_time = 0
    f = open(output_file,"w")
    for i in range(num_steps):
        f.write(f"{time_step};{num_steps};{total_time:.2f}\n")
        total_time += time_step
        time.sleep(time_step)
    f.close()


def generalized_main(fcn,argv=None,optional_arg_types={}):
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

    # Build parser
    description = dict(inspect.getmembers(fcn))["__doc__"]
    parser = argparse.ArgumentParser(prog=f"{fcn.__name__}.py",
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

            # For type == None args, parse as integer
            if param[p].default is None:
                try:
                    arg_type = optional_arg_types[p]
                except KeyError:
                    arg_type = int

            else:
                arg_type = type(param[p].default)

            parser.add_argument(f"--{p}",
                                type=arg_type,
                                default=param[p].default)

    # Parse stats
    args = parser.parse_args(argv)

    # Call function with kwargs
    fcn(**args.__dict__)

def main(argv=None):

    generalized_main(time_waster,argv)

if __name__ == "__main__":
    main()
