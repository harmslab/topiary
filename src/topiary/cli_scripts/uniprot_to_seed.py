#!/usr/bin/env python3
"""
CLI script to convert a UniProt FASTA file into a topiary seed dataframe.
"""

import topiary
from topiary.io.uniprot_to_seed import uniprot_to_seed
from topiary._private.wrap import wrap_function

import pandas as pd
import sys

def main(argv=None):
    """
    Main function for the script.
    """

    if argv is None:
        argv = sys.argv[1:]

    # Set arg types for args with None as default
    optional_arg_types = {"overwrite": bool}

    extra_args = [("--out", {"type": str, "default": "uniprot_seed.csv"}),
                  ("--overwrite", {"action": "store_true"})]

    description = \
    """
    Convert a UniProt FASTA file into a topiary seed dataframe. This script
    parses UniProt-specific headers to extract species, gene names, and
    aliases, producing a dataframe formatted for use as a seed in topiary.
    """

    # Wrap and run function
    ret, args = wrap_function(uniprot_to_seed,
                              argv=argv,
                              optional_arg_types=optional_arg_types,
                              extra_args=extra_args,
                              description=description)

    out_file = args.__dict__["out"]
    overwrite = args.__dict__["overwrite"]
    
    # Write the dataframe to a CSV file
    topiary.write_dataframe(ret, out_file, overwrite=overwrite)
    
    print(f"Seed dataframe written to {out_file}")

if __name__ == "__main__":
    main()
