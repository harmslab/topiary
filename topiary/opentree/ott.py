"""
Get OTT ids given a topiary dataframe.
"""

import topiary
from topiary._private import check
from topiary.opentree.util import species_to_ott
from .util import species_to_ott, ott_to_resolvable

import pandas as pd
import numpy as np

import re, copy

def get_df_ott(df,verbose=True,keep_anyway=False):
    """
    Return a copy of df with an ott column holding open tree of life
    names for each species. It also adds a "resolvable" column indicating
    whether the species can be resolved on the open tree of life synthetic
    tree.

    Parameters
    ----------
    df : pandas.DataFrame
        dataframe that has an ott column with Open Tree of Life taxon ids
    verbose : bool, default=True
        whether or not to print out unresolvable taxa, etc.
    keep_anyway : bool, default=False
        Do not set keep = False for species that cannot be found or resolved on
        the OTT.

    Returns
    -------
    topiary_df : pandas.DataFrame
        Copy of df with added ott, orig_species, and resolvable columns.
        ott column holds ott index for the species. orig_species holds what used
        to be in the species column. The species column is replaced by clean
        species name used by Open Tree of Life. The resolvable column indicates
        whether the species can be resolved on the synthetic tree. Rows with
        species that have no ott and/or are not resolvable have keep set to
        False. This will *not* respect the always_keep column, and will warn
        the user it is doing so. If keep_anyway is set to True, this function
        will populate the columns as described above, but will not set keep 
        to False for bad values. 
    """

    # Make sure this is a topiary dataframe
    df = check.check_topiary_dataframe(df)

    # Make sure verbose can be properly handled
    verbose = check.check_bool(verbose,"verbose")

    # Make copy of df and copy current species to original species
    local_df = df.copy()
    local_df["orig_species"] = local_df.loc[:,"species"]

    # Get ott and clean species names from opentree database
    ott_list, species_list, ott_results_dict = species_to_ott(local_df.species)
    good_ott_mask = np.array([o is not None for o in ott_list],dtype=bool)
    new_ott = [f"ott{o}" for o in np.array(ott_list)[good_ott_mask]]

    # Get species names that were not recognized
    bad_ott_mask = np.logical_not(good_ott_mask)
    unrecognized_name = np.array(species_list)[bad_ott_mask]
    unrecognized_name = list(set(unrecognized_name))
    unrecognized_name.sort()

    # Load in otts and species
    local_df["ott"] = [pd.NA for _ in ott_list]
    local_df.loc[good_ott_mask,"ott"] = new_ott
    local_df["species"] = species_list

    # Figure out which species are resolvable
    good_mask = np.array([o is not None for o in ott_list],dtype=bool)
    ott_array = np.array(ott_list)
    resolvable = np.array(ott_to_resolvable(ott_array[good_mask]),dtype=bool)
    local_df["resolvable"] = False
    local_df.loc[good_mask,"resolvable"] = np.array(resolvable)

    # Get list of taxa that could be recognized but not placed
    resolvable = local_df.loc[:,"resolvable"]
    good_not_resolved = np.logical_and(good_ott_mask,np.logical_not(resolvable))
    unresolved_taxa = np.array(species_list)[good_not_resolved]
    unresolved_taxa = list(set(unresolved_taxa))
    unresolved_taxa.sort()

    # Get bad ott or not resolvable
    bad = np.logical_or(pd.isnull(local_df["ott"]),
                        np.logical_not(local_df["resolvable"]))

    # See if always_keep species are not resolvable. These will be dropped,
    # overriding always_keep.
    try:
        always_keep = df.loc[:,"always_keep"]

        bad_drop = np.logical_and(always_keep,bad)
        bad_rows = local_df.loc[bad_drop,:]

        # Don't do this check if we are keeping regardless of whether resolvable
        if not keep_anyway:

            if len(bad_rows) > 0:
                w = "Not all entries with always_keep == True are resolvable on\n"
                w += "the opentreeoflife synthetic tree. The following rows will\n"
                w += "have always_keep and keep both set to False.\n"
                for idx in bad_rows.index:
                    uid = bad_rows.loc[idx,"uid"]
                    species = bad_rows.loc[idx,"species"]
                    ott = bad_rows.loc[idx,"ott"]
                    w += f"    {uid}: {species} ({ott})\n"
                print(w)

            local_df.loc[bad_drop,"always_keep"] = False

    except KeyError:
        pass

    # Set taxa that could not be resolved to False
    if not keep_anyway:
        local_df.loc[bad,"keep"] = False

    # Print warning data for user -- species we could not find OTT for
    if len(unrecognized_name) != 0:

        w = "\n"
        w += "Could not find OTT for following species:\n"
        for u in unrecognized_name:
            w += f"    {u}\n"

        w += "\nSetting `keep = False` for all of these species\n"
        w += re.sub("        ","",
        """
        topiary looks up unique identifiers for every species (OTT ids) using the
        opentreeoflife database. This did not work for the species listed above.
        For the moment, topiary has simply set `keep = False` for any sequences
        from these species in the dataframe, meaning they will be excluded from
        the analysis. If you want to keep these sequences, you can look up the
        OTT for the species manually on https://tree.opentreeoflife.org/,
        manually add that ott to the dataframe, and set `keep = True` for that
        row.

        This is often caused when a species has two names (for example,
        Apteryx mantelli mantelli vs. Apteryx australis mantelli). If ncbi uses
        one species name and opentreeoflife uses another, this will lead to this
        error.  This can also occur when the ncbi species name is ambiguous,
        referring to genus/species pair that has more than one subspecies
        annotated in opentreeoflife. A final common problem is when the ncbi
        sequence comes from a hybrid (for example, Bos indicus x Bos taurus).
        This is a unique species, but can't be placed on a bifurcating species
        tree.

        If you are able to find a name for the species that successfully resolves
        on the opentreeoflife database, you can update the dataframe. For the
        example of Apteryx mantelli mantellii above, you could fix this error
        by running the following code. (Note we set `keep = True` because the
        failed look up automatically set `keep = False` fo these species.)

        ```
        df.loc[df.species == "Apteryx mantelli mantelli","species"] = "Apteryx australis mantelli"
        df.loc[df.species == "Apteryx australis mantelli","keep"] = True
        df = topiary.get_df_ott(df)
        ```
        \n""")

        if verbose:
            print(w)

    # Print warning data for user -- species we could not resolve
    if len(unresolved_taxa) != 0:

        w = "\n"
        w += "Following species have OTT, but cannot be placed on tree:\n"
        for u in unresolved_taxa:
            w += f"    {u}\n"
        w += "\nSetting `keep = False` for all of these species.\n"
        w += re.sub("        ","",
        """
        topiary looks up unique identifiers for every species (OTT ids) using the
        opentreeoflife database. This did not work for the species listed above.
        For the moment, topiary has simply set `keep = False` for any sequences
        from these species in the dataframe, meaning they will be excluded from
        the analysis.

        This particular problem usually occurs for hybrids (i.e. mules). These
        cannot be placed on a branching tree because they are the result of a
        cross between branches. You have two options to fix this problem. 1)
        Remove the sequence from the analysis. This is what will happen if you
        leave `keep = False`. 2) If you need this sequence, look up the OTT for
        one of the hybrid parents (e.g. horse/donkey for mule) on
        https://tree.opentreeoflife.org/ and then use that as the OTT for this
        sequence. This is only recommended if you do not have the parent species
        in the alignment already. After adding the new OTT, make sure to set
        `keep = True` for this sequence.
        """)

        if verbose:
            print(w)

    return local_df
