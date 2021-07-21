from opentree import OT, taxonomy_helpers
import dendropy as dp

import pandas as pd
import numpy as np

from tqdm.auto import tqdm

import re, copy


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
    unrecognized_name = set(unrecognized_name)
    if len(unrecognized_name) != 0:

        print()
        print("Could not find OTT for following species:")
        for u in unrecognized_name:
            print(f"    {u}")

        print()
        print("Setting `keep = False` for all of these species.")
        print()

    # Print warning data for user -- species we could not resolve
    unresolved_taxa = set(unresolved_taxa)
    if len(unresolved_taxa) != 0:

        print()
        print("Following species have OTT, but cannot be placed on tree:")
        for u in unresolved_taxa:
            print(f"    {u}")

        print()
        print("Setting `keep = False` for all of these species.")
        print()

    return local_df


def get_species_tree(df):
    """
    Return a dendropy cladogram of species in tree.

    df: dataframe that has an ott column with Open Tree of Life taxon ids

    Returns a dendropy tree
    """

    # Only get keep = True and get unique ott
    local_df = df.loc[df.keep,:]

    # Get only rows with unique species names
    local_df = local_df.loc[local_df.loc[:,"ott"].drop_duplicates().index,:]

    # Make sure every species has ott
    if np.sum(pd.isnull(local_df.loc[:,"ott"])) > 0:
        err = "not all species OTT in dataframe."
        raise ValueError(err)

    ret = taxonomy_helpers.labelled_induced_synth(ott_ids = list(local_df.loc[:,"ott"]),
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
    final_tree = stripped_tree.extract_tree_with_taxa_labels(taxon_names)

    # Rename labels on tree so they are ott
    for n in final_tree.taxon_namespace:
        n.label = n.label.split("ott")[1]

    return final_tree

def build_species_corrected_gene_tree(df,species_tree,gene_tree_string):
    """
    Build a species-corrected gene tree.

    df: data frame with paralog, ott, keep, and uid columns
    species_tree: species tree as dendropy object ott as leaf labels
    gene_tree_string: string newick representation of desired gene tree.
                      for example, '(((A9,A8),A12),MRP126);'

    The paralogs in the dataframe must exactly match the paralogs shown
    in the gene tree.  The ott in the data frame must be a subset of the
    ott in the species tree.  Each paralog in the data frame can only have
    each species once, meaning that having two S100A9 for Mus musculus
    is not allowed. Not all species need to be seen for each paralog.

    Returns a dendropy tree with uid as tip labels.
    """

    # Create local data frame only including 'keep = True'
    local_df = df.loc[df.keep==True,:]

    # Make sure gene_tree_string is parsable and print back out so user can
    # see what they specified
    gene_tree = dp.Tree()
    gene_tree = gene_tree.get(data=gene_tree_string,schema="newick")
    for edge in gene_tree.postorder_edge_iter():
        edge.length = 0.1

    print("Using this paralog tree")
    print(gene_tree.as_ascii_plot(width=60))
    print()

    # Get all paralogs from the specified gene tree
    paralogs_in_gene_tree = []
    for leaf in gene_tree.leaf_node_iter():
        if leaf.taxon is None:
            err = "\nall paralogs must have names in the gene tree\n"
            raise ValueError(err)
        paralogs_in_gene_tree.append(leaf.taxon.label)

    # Make sure paralogs in the gene tree are unique
    if len(paralogs_in_gene_tree) != len(set(paralogs_in_gene_tree)):
        err = "\nparalogs in gene tree must be unique\n\n"
        raise ValueError(err)

    # Get paralogs seen in gene tree as unique set
    paralogs_in_gene_tree = set(paralogs_in_gene_tree)

    # Get set of unique paralogs seen in data frame
    try:
        paralogs_in_df = set(local_df.paralog.drop_duplicates())
    except KeyError:
        err = "data frame must have a 'paralog' column.\n"
        raise ValueError(err)

    # Make sure exactly the same set of paralogs in both gene tree and df
    if paralogs_in_gene_tree != paralogs_in_df:
        err = "\nparalogs must be identical in gene tree and df\n"
        err += f"df paralogs: {paralogs_in_df}\n"
        err += f"gene tree paralogs: {paralogs_in_gene_tree}\n"
        err += "\n\n"
        raise ValueError(err)

    # Make sure each species has its own unique paralog assigned
    duplicates = {}
    for paralog in paralogs_in_df:
        species = local_df[local_df.paralog == paralog].species
        duplicate_species = list(species[species.duplicated()])
        if len(duplicate_species) > 0:
            duplicates[paralog] = duplicate_species

    if len(duplicates) > 0:
        err = "each paralog can only be assigned once to each species\n"
        err = "duplicates:\n"
        for k in duplicates:
            for this_duplicate in duplicates[k]:
                err += f"    {k} {this_duplicate}\n"
            err += "\n"
        err += "\n\n"
        raise ValueError(err)

    # Get all ott seen in the species tree
    ott_in_species_tree = []
    for leaf in species_tree.leaf_node_iter():
        ott_in_species_tree.append(leaf.taxon.label)
    ott_in_species_tree = set(ott_in_species_tree)

    # Get all ott seen in df
    ott_in_df = set([f"{ott}" for ott in local_df.ott])

    # Make sure all ott in data frame are in species tree
    if not ott_in_df.issubset(ott_in_species_tree):
        err = "every ott in data frame must be seen in species tree\n"
        raise ValueError(err)

    # Now build newick
    final_tree_str = copy.copy(gene_tree_string)
    for paralog in paralogs_in_df:

        # Make df with only this paralog
        tmp_df = local_df[local_df.paralog == paralog]

        # Dictionary mapping the ott to uid (for this paralog)
        ott_to_uid = {}
        for i in range(len(tmp_df)):
            uid = tmp_df.iloc[i].uid
            ott = tmp_df.iloc[i].ott

            ott_to_uid[f"{ott}"] = uid

        # Figure out what species need to be removed
        paralog_ott = set([f"{ott}" for ott in tmp_df.ott])
        paralog_to_trim  = ott_in_species_tree - paralog_ott

        # Remove species from paralog tree
        paralog_tree = copy.deepcopy(species_tree)
        paralog_tree.prune_taxa_with_labels(paralog_to_trim)

        # Rename leaves on paralog tree to uid
        for leaf in paralog_tree.leaf_node_iter():
            leaf.taxon.label = ott_to_uid[leaf.taxon.label]

        # Write out paralog tree as a newick string
        paralog_tree_str = r"{}".format(paralog_tree.as_string(schema="newick")[:-2])

        # Replace the gene name with the new tree
        final_tree_str = re.sub(paralog,paralog_tree_str,final_tree_str)

    # Generate a final dendropy tree
    final_tree = dp.Tree()
    final_tree = final_tree.get(data=final_tree_str,
                                schema="newick")
    for edge in final_tree.postorder_edge_iter():
        edge.length = 0.01

    return final_tree
