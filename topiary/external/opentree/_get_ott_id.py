__description__ = \
"""
Get OTT ids given a topiary dataframe.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-07"

import topiary
from topiary import _arg_processors
from topiary.external.opentree.primitive import species_to_ott

from opentree import OT, taxonomy_helpers
import dendropy as dp
import ete3

import pandas as pd
import numpy as np

import re, copy

def get_ott_id(df,
               phylo_context="All life"):
    """
    Return a copy of df with an ott column holding open tree of life
    names for each species.

    Parameters
    ----------
        df: dataframe that has an ott column with Open Tree of Life taxon ids

        phylo_context: string. used to limit species seach for looking up species
                       ids on open tree of life.  To get latest strings recognized
                       by the database, use the following code:

                       ```
                       from opentree import OT
                       print(OT.tnrs_contexts().response_dict)
                       ```

                       As of 2021-08-16, the following are recognized. You can use
                       either the keys or values in this dictionary.

                       {'ANIMALS': ['Animals','Birds','Tetrapods','Mammals',
                                    'Amphibians','Vertebrates','Arthropods',
                                    'Molluscs','Nematodes','Platyhelminthes',
                                    'Annelids','Cnidarians','Arachnids','Insects'],
                        'FUNGI': ['Fungi', 'Basidiomycetes', 'Ascomycetes'],
                        'LIFE': ['All life'],
                        'MICROBES': ['Bacteria','SAR group','Archaea','Excavata',
                                     'Amoebozoa','Centrohelida','Haptophyta',
                                     'Apusozoa','Diatoms','Ciliates','Forams'],
                        'PLANTS': ['Land plants','Hornworts','Mosses','Liverworts',
                                   'Vascular plants','Club mosses','Ferns',
                                   'Seed plants','Flowering plants','Monocots',
                                   'Eudicots','Rosids','Asterids','Asterales',
                                   'Asteraceae','Aster','Symphyotrichum',
                                   'Campanulaceae','Lobelia']}

    Return
    ------
        Copy of df with added ott and orig_species column.  ott column
        holds ott index for the species. orig_species holds what used to be
        in the species column. The species column is replaced by the uniuqe
        species name used by Open Tree of Life.
    """

    # Make sure this is a topiary dataframe
    df = _arg_processors.process_topiary_dataframe(df)


    # Make copy of df and copy current species to original species
    local_df = df.copy()
    local_df["orig_species"] = local_df.loc[:,"species"]

    # Get unique list of species, stripping any leading/trailing spaces.
    species_list = list(local_df.species.drop_duplicates())

    # Actually get species ott from opentree database
    results, not_resolved = species_to_ott(species_list=species_list,
                                           phylo_context=phylo_context)

    # Create new, empty column for "ott" in the local df
    local_df["ott"] = pd.array(["" for _ in range(len(local_df))])

    # Go through the local_df and populate species, ott, and keep
    unrecognized_name = []
    unresolved_taxa = []
    final_ott = []
    for i in range(len(local_df)):

        # Get row name
        row_name = local_df.iloc[i].name

        # Get species and what keep status was before this move
        s = local_df.loc[row_name,"orig_species"]
        keep = local_df.loc[row_name,"keep"]

        # Default values
        ott_id = None
        species = s

        # Try to grab parsed results. If fails, no hit. Set ott_id to None
        # and keep species as is
        try:
            ott_id = results[s][0]
            species = results[s][1]
        except KeyError:
            pass

        # If ott is none, set to keep = False, and record we did not find the
        # species
        if ott_id is None:
            keep = False
            unrecognized_name.append(s)

        # If ott_id was something not resovled in the syntehci tree, set
        # keep = False and record we could not resolve
        if ott_id in not_resolved:
            keep = False
            unresolved_taxa.append(s)

        # Update the local_df keep, species, and ott
        local_df.loc[row_name,"keep"] = keep
        local_df.loc[row_name,"species"] = species

        if ott_id is None:
            local_df.loc[row_name,"ott"] = pd.NA
        else:
            local_df.loc[row_name,"ott"] = f"ott{ott_id}"

    # Print warning data for user -- species we could not find OTT for
    unrecognized_name = set(unrecognized_name)
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

        If you are able to find a name for the spieces that successfully resolves
        on the opentreeoflife database, you can update the dataframe. For the
        example of Apteryx mantelli mantellii above, you could fix this error
        by running the following code. (Note we set `keep = True` because the
        failed look up automatically set `keep = False` fo these species.)

        ```
        df.loc[df.species == "Apteryx mantelli mantelli","species"] = "Apteryx australis mantelli"
        df.loc[df.species == "Apteryx australis mantelli","keep"] = True
        df = topiary.get_ott_id(df,phylo_context="Animals")
        ```
        \n""")

        print(w)

    # Print warning data for user -- species we could not resolve
    unresolved_taxa = set(unresolved_taxa)
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

        print(w)

    return local_df
