"""
Get a species tree given a topiary dataframe.
"""

import topiary
from topiary._private import check
from .util import ott_to_species_tree

from opentree import taxonomy_helpers
import dendropy as dp
import ete3

import pandas as pd
import numpy as np

import copy

def df_to_species_tree(df,strict=False):
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
    strict : bool, default=False
        if strict, throw ValueError if a species cannot be found on opentree

    Returns
    -------
    species_tree : ete3.Tree
        An ete3 tree with branch lengths of 1, supports of 1, and only
        tip labels. Note: any polytomies are arbirarily resolved.
    dropped : list
        list of ott corresponding to dropped sequences
    """

    # Make sure this is a clean topiary dataframe
    df = check.check_topiary_dataframe(df)
    if "ott" not in df.columns:
        err = "\ndataframe must contain an ott column. This can be generated\n"
        err += "using topiary.opentree.get_df_ott\n\n"
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
    all_kept_df = df.copy()
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

    ott_as_int = [int(o[3:]) for o in ott_ids]

    final_tree, results = ott_to_species_tree(ott_list=ott_as_int)

    # Clean up the tree
    ott_seen = []
    for n in final_tree.traverse():

        # Give leaves species name as a feature
        if n.is_leaf():

            name_key = str(copy.deepcopy(n.name))

            n.add_feature("species",ott_to_df_columns["species"][name_key][0])
            n.add_feature("ott",ott_to_df_columns["ott"][name_key][0])
            n.add_feature("uid",ott_to_df_columns["uid"][name_key])

            n.name = copy.deepcopy(n.ott)

            ott_seen.append(n.ott)

    missing = results["not_resolved"]
    if len(missing) > 0:

        bad_rows = all_kept_df.loc[all_kept_df.ott.isin(missing),:]

        err = ["\n\nNot all species could be placed on the species tree!\n"]

        for idx in bad_rows.index:
            uid = bad_rows.loc[idx,"uid"]
            species = bad_rows.loc[idx,"species"]
            ott = bad_rows.loc[idx,"ott"]
            err.append(f"    {uid}: {species} ({ott})")

        err.append("\n\n")

        err = "\n".join(err)

        if strict:
            raise ValueError(err)
        else:
            print(err)

    return final_tree, missing
