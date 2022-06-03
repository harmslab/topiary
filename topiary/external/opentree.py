import topiary
from topiary import _arg_processors

from opentree import OT, taxonomy_helpers
import dendropy as dp
import ete3

import pandas as pd
import numpy as np

import re, copy

def is_allowed_phylo_context(phylo_context):
    """
    See if a phylo_context is recognizable in opentree database.

    Parameters
    ----------
        phylo_context: string holding phylognetic context

    Return
    ------
        phylo_context, validated as a string
    """

    phylo_context = str(phylo_context)

    # Make sure the phylo_context can be recognized by the by OT
    allowed_context = OT.tnrs_contexts().response_dict
    all_allowed = []
    for k in allowed_context:
        all_allowed.extend(allowed_context[k])

    if phylo_context not in all_allowed:
        err = f"\n\nphylo_context '{phylo_context}' not recognized. Should be one of:\n\n"
        for a in all_allowed:
            err += f"    {a}\n"
        err += "\n\n"
        raise ValueError(err)

    return phylo_context

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

    # Make sure the phylo_context is recognizable by OTT
    phylo_context = is_allowed_phylo_context(phylo_context)

    # Make copy of df and copy current species to original species
    local_df = df.copy()
    local_df["orig_species"] = local_df.loc[:,"species"]

    # Get unique list of species, stripping any leading/trailing spaces.
    species_list = list(local_df.species.drop_duplicates())
    species_list = [s.strip() for s in species_list]

    # Do fuzzy match for species names
    w = OT.tnrs_match(species_list,
                      context_name=phylo_context,
                      do_approximate_matching=True)

    # Go through hits
    results = {}
    for match in w.response_dict["results"]:

        # No match
        if len(match["matches"]) == 0:
            results[match['name']] = (None,match['name'])

        # Match
        else:

            matches = match["matches"]

            # Multiple match -- print warning
            if len(matches) > 1:

                # If we have multiple fuzzy matches, where the first is not exactly
                # equal to the query, warn the user
                if matches[0]['matched_name'] != matches[0]['taxon']['unique_name']:

                    w = f"\n\nSpecies {matches[0]['matched_name']} had multiple hits. Taking first.\n"
                    w += "Matches:\n"
                    for i in range(len(matches)):
                        w += f"    {matches[i]['taxon']['unique_name']}\n"
                    print(w)

            # Record data about hit
            hit = matches[0]

            matched_name = hit["matched_name"]
            otl_name = hit["taxon"]["unique_name"]
            ott_id = hit["taxon"]["ott_id"]

            results[matched_name] = (ott_id,otl_name)


    # Record non-hits
    for unmatch in w.response_dict["unmatched_names"]:
        results[unmatch] = (None,unmatch,None)

    # List of ott numbers to query as tree
    tmp_ott_list = [results[k][0] for k in results if results[k][0] is not None]

    # Check for ott that cannot be resolved into trees
    ret = taxonomy_helpers.labelled_induced_synth(ott_ids=tmp_ott_list,
                                                  label_format="name_and_id",
                                                  inc_unlabelled_mrca=False)

    # Taxa unresolvable in synthetic tree
    not_resolved = []
    for bad in ret["unknown_ids"].keys():
        ott_id = ret["unknown_ids"][bad]["ott_id"]
        not_resolved.append(ott_id)

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


def get_species_tree(df):
    """
    Return an ete3 cladogram of species in tree.

    df: dataframe that has an ott column with Open Tree of Life taxon ids

    Returns an ete3 tree with branch lengths of 1, supports of 1, and only
    tip labels. Note: any polytomies are arbirarily resolved.
    """

    # Make sure this is a clean topiary dataframe
    df = _arg_processors.process_topiary_dataframe(df)

    # Only get keep = True
    df = df.loc[df.keep,:].copy()

    # ott_to_df_columns will be a dictionary of dictionaries. The top-level
    # dictionary is keyed to different columns (nickname, uid, etc.). Each of
    # those keys maps to a dictionary mapping ott to a tuple of values in that
    # column that have this ott. These will be attached as features to the
    # leaf nodes in the final ete3 tree.
    to_get = ["nickname","species","sequence","name","ott","uid"]
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
            for k in ott_to_df_columns:
                n.add_feature(k,ott_to_df_columns[k][name_key])

            n.name = copy.deepcopy(n.species[0])

    return final_tree
