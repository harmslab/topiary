"""
Basic interface to opentree.
"""

import topiary
from topiary._private import check

from opentree import OT, taxonomy_helpers

import pandas as pd
import numpy as np

import re, copy, time

def is_allowed_phylo_context(phylo_context):
    """
    See if a phylo_context is recognizable in opentree database.

    Parameters
    ----------
    phylo_context : str
        string holding query phylognetic context

    Returns
    -------
    phylo_context : str
        validated phylo_context as a string
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

def get_phylo_context(species_list):
    """
    Get phylogenetic context (i.e. 'All life', 'Mammals') given a list of
    species names.

    Parameters
    ----------
    species_list : list
        list of binomial species names as strings (i.e. ["Homo sapiens",
        "Mus musculus"])

    Returns
    -------
    phylo_context : str
        string giving phylogenetic context
    """

    # Clean up species list
    species_list = check.check_iter(species_list,
                                    "species_list",
                                    required_value_type=str)

    # No species, return all life
    if len(species_list) == 0:
        return "All life"

    # Make sure we can get an ott for all sequences in sequence list
    ott, not_resolved = species_to_ott(species_list,phylo_context="All life")
    for k in ott:
        if ott[k][0] is None:
            err = f"\nCould not find species '{k}' in opentree database.\n"
            err += "We cannot reliably find phylogenetic context as a result.\n\n"
            raise ValueError(err)

    # Infer the phylogenetic context from this set of species
    ot_matches = OT.tnrs_infer_context(species_list)

    # Return context
    return ot_matches.response_dict["context_name"]

def get_resolvable(ott_list,chunk_size=100):
    """
    Get whether or not taxa are resolvable on the syntehtic ott tree.

    Parameters
    ----------
    ott_list : list-like
        list of ott as integers
    chunk_size : int, default=100
        break query into N chunk_size chunks to avoid slamming ott server with
        request for a huge synthetic tree

    Returns
    -------
    resolved : list
        list of ott that can be resolved
    not_resolved : list
        list of ott that cannot be resolved
    """

    ott_list = check.check_iter(ott_list,"ott_list")
    chunk_size = check.check_int(chunk_size,"chunk_size",minimum_allowed=1)

    # Check type of ott list
    try:
        ott_list = [int(o) for o in ott_list]
    except (ValueError,TypeError):
        err = "\nott_list should be a list of integer ott values\n\n"
        raise ValueError(err)

    # Make unique
    ott_list = list(set(ott_list))

    # Nothing passed in, don't do anything else
    if len(ott_list) == 0:
        return [], []

    not_resolved = []
    ott_broken_up = [ott_list[i:i + chunk_size]
                     for i in range(0, len(ott_list), chunk_size)]
    for query in ott_broken_up:

        try:
            # Get tree
            ret = taxonomy_helpers.labelled_induced_synth(ott_ids=query,
                                                          label_format="name_and_id",
                                                          inc_unlabelled_mrca=False)
            # Taxa unresolvable in synthetic tree
            for bad in ret["unknown_ids"].keys():
                ott_id = ret["unknown_ids"][bad]["ott_id"]
                not_resolved.append(ott_id)

        except ValueError:
            # opentree throws a ValueError if none of the ott can be placed on
            # on the tree. They're all bad.
            not_resolved.extend(query)

        # Pause briefly. No DOS attack here. :)
        time.sleep(0.1)

    resolved = list(set(ott_list) - set(not_resolved))

    return resolved, not_resolved



def species_to_ott(species_list,phylo_context="All life"):
    """
    Get the OTT values for a list of species.

    Parameters
    ----------
    species_list : list
        list of binomial species names as strings (i.e. ["Homo sapiens",
        "Mus musculus"])
    phylo_context : str, default="All life"
        used to limit species seach for looking up species ids on open tree of
        life

    Returns
    -------
    results : list
        list of ott found from this search
    not_resolved : list
        list of ott that could not be resolved on the species tree even though
        they have ott.
    """

    # Check input arguments
    species_list = check.check_iter(species_list,
                                    "species_list",
                                    is_not_type=str,
                                    required_value_type=str)

    # Make sure the phylo_context is recognizable by OTT
    phylo_context = is_allowed_phylo_context(phylo_context)

    # Strip leading/trailing spaces
    species_list = [s.strip() for s in species_list]

    # Do fuzzy match for species names
    ot_matches = OT.tnrs_match(species_list,
                               context_name=phylo_context,
                               do_approximate_matching=True)

    # Compile pattern to remove gobblygook like "(in domain bacteria)" that
    # comes down with tnrs match
    domain_pattern = re.compile(" \(.*?in domain.*?\)")

    # Go through hits
    results = {}
    for match in ot_matches.response_dict["results"]:

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

            matched_name = domain_pattern.sub("",hit["matched_name"])
            otl_name = domain_pattern.sub("",hit["taxon"]["unique_name"])
            ott_id = hit["taxon"]["ott_id"]

            results[matched_name] = (ott_id,otl_name)


    # Record non-hits
    for unmatch in ot_matches.response_dict["unmatched_names"]:
        results[unmatch] = (None,unmatch,None)

    # List of ott numbers to query as tree
    tmp_ott_list = [results[k][0] for k in results if results[k][0] is not None]

    # Don't try to get synthetic tree for empty sequence set
    if len(tmp_ott_list) == 0:
        return results, []

    # Get taxa that cannot be resolved
    resolved, not_resolved = get_resolvable(tmp_ott_list)

    return results, not_resolved
