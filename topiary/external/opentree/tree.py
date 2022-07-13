"""
Get a species tree given a topiary dataframe.
"""

import topiary
from topiary._private import check

from opentree import taxonomy_helpers
import dendropy as dp
import ete3

import pandas as pd
import numpy as np

import re, copy

def get_species_tree(df):
    """
    Return an ete3 cladogram of species in tree. The leaves on the tree will
    have the following features:

    + :code:`leaf.name`: ott as string
    + :code:`leaf.ott`: ott as string
    + :code:`leaf.species`: bionomial species name as string
    + :code:`leaf.uid`: list of all uid that have this species

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe that has an ott column with Open Tree of Life taxon
        ids

    Returns
    -------
    species_tree : ete3.Tree
        An ete3 tree with branch lengths of 1, supports of 1, and only
        tip labels. Note: any polytomies are arbirarily resolved.
    """

    # Make sure this is a clean topiary dataframe
    df = check.check_topiary_dataframe(df)
    if "ott" not in df.columns:
        err = "\ndataframe must contain an ott column. This can be generated\n"
        err += "using topiary.opentree.get_ott_id\n\n"
        raise ValueError(err)

    # Only get keep = True
    df = df.loc[df.keep,:].copy()

    # ott_to_df_columns will be a dictionary of dictionaries. The top-level
    # dictionary is keyed to different columns. Each of
    # those keys maps to a dictionary mapping ott to a tuple of values in that
    # column that have this ott. These will be attached as features to the
    # leaf nodes in the final ete3 tree.
    to_get = ["species","ott","uid"]
    ott_to_df_columns = dict([(k,{}) for k in to_get])
    for k in to_get:

        # If the dataframe doesn't have the relevant column, do not grab
        if k not in df.columns:
            ott_to_df_columns.pop(k)
            continue

        # Load values from the dataframe into the dictionary
        for i in range(len(df)):
            row = df.iloc[i]
            try:
                ott_to_df_columns[k][row.ott].append(row[k])
            except KeyError:
                ott_to_df_columns[k][row.ott] = [row[k]]

        # Convert to tuple for each of these, as ete3 requires features are
        # hashable types.
        for o in ott_to_df_columns[k]:
            ott_to_df_columns[k][o] = tuple(ott_to_df_columns[k][o])

    # Get only rows with unique ott
    df = df.loc[df.loc[:,"ott"].drop_duplicates().index,:]

    # Make sure every species has an ott
    if np.sum(pd.isnull(df.loc[:,"ott"])) > 0:
        err = "\nNot all species have ott in the dataframe.\n"
        raise ValueError(err)

    ott_ids = list(df.loc[:,"ott"])
    species = list(df.loc[:,"species"])
    ott_species_dict = {}
    for i in range(len(ott_ids)):
        ott_species_dict[ott_ids[i]] = species[i]

    ott_as_int = [o[3:] for o in ott_ids]
    ret = taxonomy_helpers.labelled_induced_synth(ott_ids=ott_as_int,
                                                  label_format="name_and_id",
                                                  inc_unlabelled_mrca=False)

    # Write out without all the ancestor junk returned by opentree
    t = ret["labelled_tree"].as_string(schema="newick",
                                       suppress_leaf_taxon_labels=False,
                                       suppress_leaf_node_labels=True,
                                       suppress_internal_taxon_labels=True,
                                       suppress_internal_node_labels=True,
                                       suppress_edge_lengths=True,
                                       suppress_rooting=False,
                                       suppress_annotations=True,
                                       suppress_item_comments=True)

    # Read in new tree
    stripped_tree = dp.Tree.get(data=t,schema="newick")
    taxon_names = [s.label for s in stripped_tree.taxon_namespace]

    # Extract tree with labels.  This will remove all of the empty singleton
    # entries that otl brings down.
    clean_tree = stripped_tree.extract_tree_with_taxa_labels(taxon_names)

    # Rename labels on tree so they are ott
    for n in clean_tree.taxon_namespace:
        label = n.label.split("ott")[-1]
        n.label = f"ott{label}"

    # Convert tree to ete3 tree.
    final_tree = ete3.Tree(clean_tree.as_string(schema="newick",
                                                suppress_rooting=True),
                           format=9)

    # Arbitrarily resolve any polytomies
    final_tree.resolve_polytomy()

    # Give every node a support of 1 and a branch length of 1
    for n in final_tree.traverse():
        if n.dist != 1:
            n.dist = 1
        if n.support != 1:
            n.support = 1

        # Give leaves species name as a feature
        if n.is_leaf():

            name_key = str(copy.deepcopy(n.name))

            n.add_feature("species",ott_to_df_columns["species"][name_key][0])
            n.add_feature("ott",ott_to_df_columns["ott"][name_key][0])
            n.add_feature("uid",ott_to_df_columns["uid"][name_key])

            n.name = copy.deepcopy(n.ott)

    return final_tree
