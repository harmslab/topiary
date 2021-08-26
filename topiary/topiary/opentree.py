from opentree import OT, taxonomy_helpers
import dendropy as dp
import ete3

import pandas as pd
import numpy as np

import re, copy, warnings

def get_ott_id(df,
               context_name="All life",
               strip_regex=["MRCA of taxa in",
                            "_species in domain .*?_",
                            "\(species in domain .*?\)"]):
    """
    Return a copy of df with an ott column holding open tree of life
    names for each species.

    df: dataframe that has an ott column with Open Tree of Life taxon ids
    context_name: limit search to specific groups ("Animals" etc.)
    strip_regex: regular expressions that will be found and replaced with ""
                 in taxon names coming off open tree of life.

    Returns copy of df with added ott and orig_species column.  ott column
    holds ott index for the species. orig_species holds what used to be
    in the species column. The species column is replaced by the uniuqe
    species name used by Open Tree of Life.
    """

    # Make sure the context_name can be recognized by the by OT
    allowed_context = OT.tnrs_contexts().response_dict
    all_allowed = []
    for k in allowed_context:
        all_allowed.extend(allowed_context[k])
    if context_name not in all_allowed:
        err = f"\n\ncontext_name '{context_name}' not recognized. Should be one of:\n\n"
        for a in all_allowed:
            err += f"    {a}\n"
        err += "\n\n"
        raise ValueError(err)

    # Compile regex
    strip_regex = [re.compile(p) for p in strip_regex]

    # Make copy of df and copy current species to original species
    local_df = df.copy()
    local_df["orig_species"] = local_df.loc[:,"species"]

    # Get unique list of species
    species_list = list(local_df.species.drop_duplicates())

    # Do fuzzy match for species names
    w = OT.tnrs_match(species_list,
                      context_name=context_name,
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

                    print(f"Species {matches[0]['matched_name']} had multiple hits. Taking first.")
                    for i in range(len(matches)):
                        print(f"    {matches[i]['taxon']['unique_name']}")
                    print()

            # Record data about hit
            hit = matches[0]

            matched_name = hit["matched_name"]
            otl_name = hit["taxon"]["unique_name"]
            for s in strip_regex:
                otl_name = s.sub("",otl_name)
                otl_name = otl_name.strip()

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
    local_df["ott"] = pd.array([None for _ in range(len(local_df))],
                               dtype=pd.Int64Dtype())

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
        local_df.loc[row_name,"ott"] = ott_id

    # Print warning data for user -- species we could not find OTT for
    ott_error_found = False
    unrecognized_name = set(unrecognized_name)
    if len(unrecognized_name) != 0:

        print()
        print("Could not find OTT for following species:")
        for u in unrecognized_name:
            print(f"    {u}")

        print("\nSetting `keep = False` for all of these species\n")

        ott_error_found = True

    # Print warning data for user -- species we could not resolve
    unresolved_taxa = set(unresolved_taxa)
    if len(unresolved_taxa) != 0:

        print()
        print("Following species have OTT, but cannot be placed on tree:")
        for u in unresolved_taxa:
            print(f"    {u}")
        print("\nSetting `keep = False` for all of these species.\n")

        ott_error_found = True

    if ott_error_found:

        print(re.sub("        ","",
        """
        topiary looks up unique identifiers for every species (OTT ids) using the
        opentreeoflife database. This did not work for the species listed above.
        For the moment, topiary has simply set `keep = False` for any sequences
        from these species in the dataframe, meaning they will be excluded from
        the analysis. If you want to keep these sequences, you can look up the
        OTT for the species manually on https://tree.opentreeoflife.org/.

        This is often caused when a species has two names (for example,
        Apteryx mantelli mantelli vs. Apteryx australis mantelli). If ncbi uses
        one species name and opentreeoflife uses another, this will lead to this
        error.  This can also occur when the ncbi species name is ambiguous,
        referring to genus/species pair that has more than one subspecies
        annotated in opentreeoflife. A final common problem is when the ncbi
        sequence comes from a hybrid (for example, Bos indicus x Bos taurus).
        This is a unique species, but can't be placed on a bifurcating species
        tree.

        If you are able to find a name for the speices that successfully resolves
        on the opentreeoflife database, you can update the dataframe. For the
        example of Apteryx mantelli mantellii above, you could fix this error
        by running the following code. (Note we set `keep = True` because the
        failed look up automatically set `keep = False` fo these species.)

        ```
        df.loc[df.species == "Apteryx mantelli mantelli","species"] = "Apteryx australis mantelli"
        df.loc[df.species == "Apteryx australis mantelli","keep"] = True
        df = topiary.get_ott_id(df,context_name="Animals")
        ```
        """))

    return local_df


def get_species_tree(df):
    """
    Return an ete3 cladogram of species in tree.

    df: dataframe that has an ott column with Open Tree of Life taxon ids

    Returns an ete3 tree with branch lengths of 1, supports of 1, and only
    tip labels. Note: any polytomies are arbirarily resolved.
    """

    # Only get keep = True
    df = df.loc[df.keep,:].copy()

    # coerce ott to string int
    tmp_ott = [f"{o}" for o in np.array(df.loc[:,"ott"],dtype=np.int64)]
    df.loc[:,"ott"] = tmp_ott

    # Get only rows with unique ott
    df = df.loc[df.loc[:,"ott"].drop_duplicates().index,:]

    # Make sure every species has ott
    if np.sum(pd.isnull(df.loc[:,"ott"])) > 0:
        err = "not all species OTT in dataframe."
        raise ValueError(err)

    ott_ids = list(df.loc[:,"ott"])
    species = list(df.loc[:,"species"])
    ott_species_dict = {}
    for i in range(len(ott_ids)):
        ott_species_dict[ott_ids[i]] = species[i]

    ret = taxonomy_helpers.labelled_induced_synth(ott_ids = ott_ids,
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
        n.label = label

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
            n.add_feature("species",ott_species_dict[n.name])

    return final_tree
