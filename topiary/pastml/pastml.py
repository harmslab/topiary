"""
Interface to pastml library.
"""

import pastml
from pastml import acr

import ete3

import numpy as np
import pandas as pd

def get_ancestral_gaps(alignment_file,tree_file,prediction_method="DOWNPASS"):
    """
    Get ancestral gaps from an alignment and raxml output tree file. Gaps are
    reconstructed by parsimony using the DOWNPASS algorithm as implemented in
    pastml.

    Parameters
    ----------
    alignment_file : str
        phy file used to generate ancestors in RAxML
    tree_file : str
        output tree file with labeled internal nodes
    prediction_method : str, default="DOWNPASS"
        method to reconstruct gaps.

    Returns
    -------
    gap_anc_dict : dict
        dictionary keying internal node names to lists of True (gap), False (no
        gap), and None (gapping unclear) for each site in that ancestor.
    """

    # Read the alignment file
    counter = 0
    leaf_names = []
    with open(alignment_file) as f:
        for line in f:

            # First line
            if counter == 0:
                col = line.split()
                num_taxa = int(col[0])
                num_sites = int(col[1])
                counter += 1

                char_matrix = np.zeros((num_taxa,num_sites),dtype=np.uint8)
                continue

            # Next, blank line
            if line.strip() == "":
                counter += 1
                continue

            # Alternating leaf id and sequence lines
            if counter % 2 == 0:
                leaf_names.append(line.strip())
                counter += 1
            else:
                index = (counter - 3)//2
                char_matrix[index,:] = np.array([c == "-" for c in line.strip()])
                counter += 1

    # Create a data frame where indexes are leaf names and columns are each gap
    out_dict = {}
    for i in range(char_matrix.shape[1]):
        out_dict[f"g{i}"] = char_matrix[:,i]
    gap_df = pd.DataFrame(out_dict)
    gap_df.index = leaf_names

    # Gaps, named as column names
    gap_names = list(gap_df.columns)

    # Load the tree, keeping the internal node names
    tree = ete3.Tree(tree_file,format=1)

    # Reconstruct gaps across tree
    acr_result = acr.acr(tree,gap_df,prediction_method=prediction_method)

    # Create dictionary keying ancestor name to gap status across sequence
    gap_anc_dict = {}
    for node in tree.traverse('preorder'):
        if node.name in leaf_names:
            continue

        gap_anc_dict[node.name] = []

        # Go through gaps and grab from node feature
        for g in gap_names:
            state = node.__dict__[g]
            if len(state) == 1:
                if state == {0}:
                    state = False
                else:
                    state = True
            else:
                state = None

            gap_anc_dict[node.name].append(state)

    return gap_anc_dict
