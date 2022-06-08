__description__ = \
"""
Basic interface to opentree.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-07"

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

def get_phylo_context(species_list):
    """
    Get phylogenetic context (i.e. 'All life', 'Mammals') given a list of
    species names.

    Parameters
    ----------
        species_list: list of species as strings

    Return
    ------
        string giving phylogenetic context
    """

    # Clean up species list
    species_list = _arg_processors.process_iter(species_list,
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
    w = OT.tnrs_infer_context(species_list)

    # Return context
    return w.response_dict["context_name"]

def species_to_ott(species_list,phylo_context="All life"):
    """
    Get the OTT values for a list of species.

    Parameters
    ----------
        species_list: list of species as strings
        phylo_context: string. used to limit species seach for looking up species
                       ids on open tree of life

    Return
    ------
        ott results for search, list of ott that could not be resolved even
        though they have ott.
    """

    # Check input arguments
    species_list = _arg_processors.process_iter(species_list,
                                                "species_list",
                                                required_value_type=str)

    # Make sure the phylo_context is recognizable by OTT
    phylo_context = is_allowed_phylo_context(phylo_context)

    # Strip leading/trailing spaces
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

    # Don't try to get synthetic tree for empty sequence set
    if len(tmp_ott_list) == 0:
        return results, []

    # Check for ott that cannot be resolved into trees
    ret = taxonomy_helpers.labelled_induced_synth(ott_ids=tmp_ott_list,
                                                  label_format="name_and_id",
                                                  inc_unlabelled_mrca=False)

    # Taxa unresolvable in synthetic tree
    not_resolved = []
    for bad in ret["unknown_ids"].keys():
        ott_id = ret["unknown_ids"][bad]["ott_id"]
        not_resolved.append(ott_id)

    return results, not_resolved
