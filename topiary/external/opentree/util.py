"""
Functions to interact directly with opentree database
"""

import topiary
from topiary._private import check

import opentree
from opentree import OT, taxonomy_helpers
import dendropy as dp
import ete3

import pandas as pd
import numpy as np

import re, copy, time

def _validate_ott_vs_species(ott_list=None,species_list=None):
    """
    Take the ott_list and species_list input, validate, and return an ott_list.

    Parameters
    ----------
    ott_list : list or None
        input ott_list to check
    species_list : list or None
        input species_list to check

    Returns
    -------
    ott_list : list
        validated ott list
    """

    failed_combo = False

    if ott_list is None and species_list is None:
        failed_combo = True

    if ott_list is not None and species_list is not None:
        failed_combo = True

    if failed_combo:
        err = "\nEither ott_list or species_list must be specified, but not both.\n\n"
        raise ValueError(err)

    if species_list is not None:
        ott_list = species_to_ott(species_list)[0]

    ott_list = check.check_iter(ott_list,"ott_list",required_value_type=int)

    return ott_list


def species_to_ott(species):
    """
    Return ott ids (and other information) given a list of species.

    Parameters
    ----------
    species : list
        list of species in binomial format

    Returns
    -------
    ott_list : list
        list of ott as integer for all species in the order they were passed.
        If species not found, set to None.
    species_list : list
        list of species found for each species. If not found, return the input
        species name.
    results : dict
        dictionary of information about species keyed to input species name

    Notes
    -----
    For all species, whether matched or not, the results dictionary will have
    following fields keyed to input species name.

     + matched: whether or not this gave an unambiguous match
     + num_matches: number of matched hits
     + msg: message describing information about match
     + ret: match tuple from opentree.OT.tnrs_match
     + ott_id: integer ott id
     + ott_name: found name corresponding to ott_id
     + taxid: NCBI taxid or None if no NCBI taxid associated with record
     + resolved: bool. whether or not species is resolved on synthetic tree

    The opentree tnrs_match api is smart enough to infer context from
    the species you pass in. If you pass in an ambiguous species name, it
    will select the species based on the other species in the list. If the
    species cannot be inferred from the other species context, this function
    will return no ott for the ambiguous species.
    """

    # Make sure species list is sane
    species = topiary._private.check.check_iter(species,
                                                "species",
                                                required_value_type=str,
                                                is_not_type=str)

    # Get unique species.
    unique_species = list(set(species))
    unique_species = [s.strip() for s in species]

    # Grab species names
    ot_matches = OT.tnrs_match(unique_species,do_approximate_matching=True)

    results = {}
    for i, result in enumerate(ot_matches.response_dict["results"]):

        s = unique_species[i]

        # Parse matches: 0, 1, or many
        matches = result["matches"]
        if len(matches) == 0:
            msg = f"Not match for '{s}'.\n"
            results[s] = {"matched":False,
                          "num_matches":0,
                          "msg":msg,
                          "ret":matches,
                          "ott_id":None,
                          "ott_name":None,
                          "taxid":None}
            continue

        elif len(matches) == 1:
            taxon = matches[0]["taxon"]
            results[s] = {"matched":not matches[0]["is_approximate_match"],
                          "num_matches":1,
                          "msg":"success",
                          "ret":matches}
        else:
            taxon = matches[0]["taxon"]

            msg = "success"
            if matches[0]["is_approximate_match"]:
                msg = f"No exact match for '{s}'. Approximate matches:\n"

                for m in matches:
                    msg += f"    {m['matched_name']}\n"
                msg += "\n"

            results[s] = {"matched":not matches[0]["is_approximate_match"],
                          "num_matches":len(matches),
                          "msg":msg,
                          "ret":matches}

        results[s]["ott_id"] = taxon["ott_id"]
        results[s]["ott_name"] = taxon["name"]

        # Try to get taxid from the ott taxon entry
        taxid = None
        try:
            tax_sources = taxon["tax_sources"]
            for t in tax_sources:
                if t.startswith("ncbi"):
                    taxid = int(t.split(":")[1])
                    break
        except KeyError:
            pass

        if taxid is None:
            try:
                taxid = topiary.ncbi.entrez.get_taxid(s)
            except RuntimeError:
                pass

        results[s]["taxid"] = taxid

    # Create list of ott from these results
    ott_list = []
    species_list = []
    for s in species:
        m = results[s]
        if m["msg"] == "success":
            ott_list.append(m["ott_id"])
            species_list.append(m["ott_name"])
        else:
            ott_list.append(None)
            species_list.append(s)

    return ott_list, species_list, results

def ott_species_tree(ott_list=None,species_list=None):
    """
    Get a species tree from a list of ott.

    Parameters
    ----------
    ott_list : list, optional
        list of ott ids (integers). this or species_list must be specified
    species_list : list, optional
        list of binomial species (str). this or ott_list must be specified

    Returns
    -------
    species_tree : ete3.Tree or None
        species tree. None if tree cannot be pulled down.
    results : dict
        dictionary with resolved, missing, not_resolved, and not_monophyletic
        ott.
    """

    ott_list = _validate_ott_vs_species(ott_list,species_list)

    # Check type of ott list
    try:
        ott_list = [int(o) for o in ott_list]
    except (ValueError,TypeError):
        err = "\nott_list should be a list of integer ott values\n\n"
        raise ValueError(err)

    # If length of ott_list is zero
    if len(ott_list) == 0:

        results = {"resolved":[],
                   "not_resolved":[],
                   "unknown_ids":[],
                   "not_monophyletic":[]}

        return None, results

    # Make unique
    ott_list = list(set(ott_list))

    # Pull down the synthetic tree from the ott server
    try:
        ret = taxonomy_helpers.labelled_induced_synth(ott_ids=ott_list,
                                                      label_format="name_and_id",
                                                      inc_unlabelled_mrca=False)

    # Value error if *all* otts are bad.
    except ValueError:
        results = {"resolved":[],
                   "not_resolved":ott_list[:],
                   "unknown_ids":ott_list[:],
                   "not_monophyletic":[]}
        return None, results


    # Taxa unresolvable in synthetic tree
    unknown_ids = []
    for bad in ret["unknown_ids"].keys():
        unknown_ids.append(int(bad[3:]))

    # Non-monophyletic taxa
    not_monophyletic = []
    for bad in ret["non-monophyletic_taxa"].keys():
        ott_id = ret["non-monophyletic_taxa"][bad]["ott_id"]
        not_monophyletic.append(ott_id)

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

    # Read in the tree
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
                                                suppress_rooting=True),format=9)

    # Arbitrarily resolve any polytomies
    final_tree.resolve_polytomy()

    # Give every node a support of 1 and a branch length of 1
    ott_seen = []
    for n in final_tree.traverse():
        if n.dist != 1:
            n.dist = 1
        if n.support != 1:
            n.support = 1

        if n.is_leaf():
            ott_seen.append(int(n.name[3:]))

    # Not seen in tree
    not_resolved = list(set(ott_list)-set(ott_seen))

    results = {"resolved":ott_seen,
               "not_resolved":not_resolved,
               "unknown_ids":unknown_ids,
               "not_monophyletic":not_monophyletic}

    return final_tree, results


def ott_resolvable(ott_list=None,species_list=None):
    """
    Get whether or not taxa are resolvable on the synthetic ott tree.

    Parameters
    ----------
    ott_list : list, optional
        list of ott ids (integers). this or species_list must be specified
    species_list : list, optional
        list of binomial species (str). this or ott_list must be specified

    Returns
    -------
    resolvable : list
        list of True/False for each ott in ott_list
    """

    # Check ott list
    ott_list = _validate_ott_vs_species(ott_list,species_list)

    try:
        ott_list = [int(o) for o in ott_list]
    except (ValueError,TypeError):
        err = "\nott_list should be a list of integer ott values\n\n"
        raise ValueError(err)

    # Return empty list if empty input
    if len(ott_list) == 0:
        return []

    # Get tree and metadata
    final_tree, results = ott_species_tree(ott_list)

    # Dictionary of resolved/not resolved keyed to ott
    resolved_dict = dict([(r,True) for r in results["resolved"]])
    for r in results["not_resolved"]:
        resolved_dict[r] = False

    # Convert to a final True/False list
    resolvable = []
    for o in ott_list:
        resolvable.append(resolved_dict[o])

    return resolvable


def ott_mrca(ott_list=None,species_list=None,move_up_by=0,avoid_all_life=True):
    """
    Get the most recent common ancestor given a list of ott. Unrecognized ott
    are dropped with a warning.

    Parameters
    ----------
    ott_list : list, optional
        list of ott ids (integers). this or species_list must be specified
    species_list : list, optional
        list of binomial species (str). this or ott_list must be specified
    move_up_by : int, default=0
        starting at actual MRCA, move up by this number of ranks. For
        example, if the MRCA for a set of OTT was a kingdom and
        move_up_by = 1, this would yield the relevant domain.
    avoid_all_life : bool, default=True
        if possible, avoid the jump to all cellular organisms. This takes
        precedence over move_up_by.

    Returns
    -------
    out : dict
        dictionary with keys ott_name, ott_id, ott_rank, lineage, and taxid.
    """

    ott_list = _validate_ott_vs_species(ott_list,species_list)

    avoid_all_life = check.check_bool(avoid_all_life,
                                      "avoid_all_life")
    move_up_by = check.check_int(move_up_by,
                                 "move_up_by",
                                 minimum_allowed=0)

    if len(ott_list) == 0:
        out = {}
        out["ott_name"] = "life"
        out["ott_id"] = 805080
        out["ott_rank"] = "no rank"
        out["lineage"] = None
        out["taxid"] = 1

        return out

    # Get the synthetic mrca
    try:
        synth_mrca = OT.synth_mrca(ott_ids=ott_list)
    except opentree.OTWebServicesError as e:
        err = "No valid OTT in ott list.\n"
        raise ValueError(err) from e

    # Get ott for the mrca species
    mrca = synth_mrca.response_dict["mrca"]
    try:
        taxon = mrca["taxon"]
    except KeyError:
        taxon = synth_mrca.response_dict["nearest_taxon"]

    # Get ott, name, and lineage for the mrca
    mrca_ott = taxon["ott_id"]
    mcra_name = taxon["name"]
    mrca_info = OT.taxon_info(ott_id=mrca_ott,include_lineage=True)
    lineage = mrca_info.response_dict["lineage"]

    # Start lineage with the mrca
    lineage.insert(0,taxon)

    # Whack off "life" if that's the highest rank and we're avoiding all life
    if len(lineage) > 1 and avoid_all_life:
        if lineage[-1]["name"] == "life":
            lineage = lineage[:-1]

    # Whack off "cellular life" if that's the highest rank and we're avoiding
    # all life
    if len(lineage) > 1 and avoid_all_life:
        if lineage[-1]["name"] == "cellular organisms":
            lineage = lineage[:-1]

    # Figure out how much to move up by given what's in lineage
    if move_up_by > len(lineage) - 1:
        move_up_by = len(lineage) - 1

    anc = lineage[move_up_by]

    # Get information about this ancestor
    out = {}
    out["ott_name"] = anc["name"]
    out["ott_id"] = anc["ott_id"]
    out["ott_rank"] = anc["rank"]
    out["lineage"] = lineage
    taxid = None
    for t in anc["tax_sources"]:
        if t.startswith("ncbi"):
            taxid = int(t.split(":")[1])
            break

    if taxid is None:
        try:
            taxid = topiary.ncbi.entrez.get_taxid(anc["name"])
        except RuntimeError:
            pass

    out["taxid"] = taxid

    return out
