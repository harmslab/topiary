#!/usr/bin/env python3
"""
Command line interface to reports.create_report
"""

import topiary
from topiary.reports import tree_report
from topiary._private.wrap import wrap_function

import sys

def main(argv=None):

    if argv is None:
        argv = sys.argv[1:]

    # Set arg types for args with None as default
    optional_arg_types = {"ancestor_directory":str}

    # Wrap and run function
    wrap_function(tree_report,
                  argv=argv,
                  optional_arg_types=optional_arg_types)

if __name__ == "__main__":
    main()
