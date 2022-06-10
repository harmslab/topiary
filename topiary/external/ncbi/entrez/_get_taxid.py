__description__ = \
"""
Get the NCBI taxid for a list of species.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-08"

import topiary
from topiary import _arg_processors

from Bio import Entrez

def get_taxid(species_list):
    """
    Return the taxid for a list of species. NOTE: these are not guaranteed to
    be in the same order as the input species.

    Parameters
    ----------
        species_list: list of species in binomial format (i.e. Homo sapiens).

    Return
        list of taxid.
    """

    # Make sure species list is sane, each species is unique, and the list
    # is sorted
    species_list = _arg_processors.process_iter(species_list,
                                                "species_list",
                                                required_value_type=str)
    return_singleton = False
    if type(species_list) is str:
        species_list = [species_list]
        return_singleton = True

    species_list = list(set(species_list))
    species_list.sort()

    # return nothing if nothing was passed in
    if len(species_list) == 0:
        return []

    # Create a search term
    search_term = " OR ".join(species_list)

    # Access Entrez and download
    handle = Entrez.esearch(db="taxonomy",
                            retmax=len(species_list)*2,
                            term=search_term,
                            idtype="uilist")
    record = Entrez.read(handle)

    # See how many records we pulled down
    count = int(record["Count"])

    if count != len(species_list):

        error_list = []

        error = record["ErrorList"]
        for k in error:
            if len(error[k]) > 0:
                error_list.append(f"{k}: {error[k]}")

        err = "\nEntrez returned the following error when retrieving taxids. This\n"
        err += "caused by mis-spelled species name(s).\n\n"
        for e in error_list:
            err += f" -> {e}\n"
        err += "\n"

        raise RuntimeError(err)

    taxid_list = [str(i) for i in record["IdList"]]

    if return_singleton:
        taxid_list = taxid_list[0]

    return taxid_list
