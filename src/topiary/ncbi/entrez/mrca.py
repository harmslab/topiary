

"""
Use entrez to get the MRCA taxid for a list of species.
"""

import topiary
from topiary._private import check

from Bio import Entrez
import re

def get_mrca_taxid(species_list):
    """
    Use entrez to get the MRCA taxid for a list of species.

    Parameters
    ----------
    species_list : list
        list of species in binomial format (i.e. Homo sapiens).

    Returns
    -------
    mrca_taxid : int
        NCBI taxid of the most recent common ancestor.
    """

    # Get taxids for all species
    taxids = topiary.ncbi.get_taxid(species_list)

    if len(taxids) == 0:
        return 1 # Root of all life

    if len(taxids) == 1:
        # Get the lineage for this single taxid to find its parent or just return it?
        # Usually MRCA of one thing is itself.
        return int(taxids[0])

    # Fetch lineages for all taxids
    handle = Entrez.efetch(db="taxonomy", id=",".join(taxids), retmode="xml")
    records = Entrez.read(handle)

    lineages = []
    for record in records:
        lineage = [1] # Start with root
        for taxon in record["LineageEx"]:
            lineage.append(int(taxon["TaxId"]))
        lineage.append(int(record["TaxId"]))
        lineages.append(lineage)

    if not lineages:
        return 1

    # Find the last common taxid in all lineages
    common_ancestor = 1
    for i in range(min(len(l) for l in lineages)):
        current_taxid = lineages[0][i]
        if all(l[i] == current_taxid for l in lineages):
            common_ancestor = current_taxid
        else:
            break

    return common_ancestor
