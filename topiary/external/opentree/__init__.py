__description__ = \
"""
Interface to open tree of life database.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-07"

from .primitive import is_allowed_phylo_context, get_phylo_context, species_to_ott
from ._get_ott_id import get_ott_id
from ._get_species_tree import get_species_tree

import pandas as pd
import os

# -----------------------------------------------------------------------------
# Load a csv file mapping between open tree phylogenetic context and ncbi
# taxid

def _load_phylo_to_taxid(csv):
    """
    Load the mapping between open tree of life phylogenetic contexts and NCBI
    taxids.

    Parameters
    ----------
        csv: csv file with mapping

    Return
    ------
        dictionary with mapping
    """

    # Read csv
    df = pd.read_csv(csv)

    # Go through dataframe
    phylo_to_taxid = {}
    for idx in df.index:

        # Cast context and taxid as strings
        context = str(df.loc[idx,"phylo_context"])
        taxid = str(df.loc[idx,"taxid"]).strip()

        # Pull out "None" taxid
        if taxid == "None":
            taxid = None
            phylo_to_taxid[context] = taxid
            continue

        # Pull out tuple taxid
        if taxid.startswith("("):
            chopped = taxid[1:-1].split(",")
            taxid = tuple([int(c) for c in chopped])
            phylo_to_taxid[context] = taxid
            continue

        # Grab simple integer taxid
        taxid = int(taxid)
        phylo_to_taxid[context] = taxid

    # Update keys to include uppercase and lowercase versions
    keys = list(phylo_to_taxid.keys())
    for k in keys:
        phylo_to_taxid[k.lower()] = phylo_to_taxid[k]
        phylo_to_taxid[k.upper()] = phylo_to_taxid[k]

    return phylo_to_taxid

# Get location of topiary in the file system
_location = os.path.dirname(os.path.realpath(__file__))
_context_to_taxid_file = os.path.join(_location,"phylo_context_to_taxid.csv")

# Try to read csv file with the ott context to taxid map
try:
    phylo_to_taxid = _load_phylo_to_taxid(_context_to_taxid_file)
except Exception as e:
    err = f"\nCould not load {_context_to_taxid_file}\n\n"
    raise Exception from e(err)
