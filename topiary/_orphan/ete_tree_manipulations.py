__description__ = \
"""
Some handy functions for manipulating ete3 trees. I developed as part of a
pipeline, then changed strategies and so these were orphaned. Good code; likely
helpful in future.
"""
__date__ = "2021-08-06"

def copy_node_features(copy_from,copy_to):
    """
    Copy features from one node to another. Use copy.deepcopy to make sure
    complex features are preserved. This copies features from copy_from
    and wipes them out on copy_to. Any features unique to copy_to are
    preserved. Note this works on the object itself, not a copy.

    copy_from: ete3 TreeNode object with features to copy
    copy_to: ete3 TreeNode object with features to load in

    returns copy_to with newly copied features
    """

    for f in copy_from.features:
        try:
            copy_to.add_feature(f,copy.deepcopy(copy_from.__dict__[f]))
        except KeyError:
            pass

    return copy_to

def local_species_correct(graft_point,clade_anc,species_tree,species_to_grab):


    species_subtree = species_tree.copy(method="deepcopy")
    species_subtree.prune(species_to_grab)

    # Copy information from left leaves onto left species tree.
    for leaf in clade_anc.get_leaves():
        match_ott = species_subtree.search_nodes(**{"name":leaf.ott})
        if len(match_ott) > 1:
            err = "species tree not unique!\n"
            raise RuntimeError(err)
        if len(match_ott) == 0:
            err = "species in gene tree but not species tree\n"
            raise RuntimeError(err)

        species_leaf = match_ott[0]
        species_leaf = copy_node_features(leaf,species_leaf)

    # Graft on species tree with information copied in from gene tree
    graft_point.add_child(species_subtree)

    # Hack off original gene tree clade
    clade_anc.detach()



def merge_sisters(some_tree,work_on_copy=True):
    """
    Merge directly sister duplicates. This would merge the top three A
    because they are directly sister, but leave the bottom A intact.
    merged nodes have the feature merge_stack populated, which holds
    a dictionary of features for the merged nodes.

           _A             ___A
         _/              /
       _/ \_A          _/   _B
      | \___A     ->    \__/
     -|    _B              \_A
      |___/
          \_A

    some_tree: ete3 tree with duplicates annotated using annotate_duplicates
    work_on_copy: work on a copy rather than the source tree

    returns tree with sister duplicates merged
    """

    # Copy tree -- do not operate on input tree directly
    if work_on_copy:
        tree = some_tree.copy(method="deepcopy")
    else:
        tree = some_tree

    # Create list of duplicate (species|paralog) calls
    duplicate_calls = {}
    for node in tree.get_leaves():
        try:
            duplicate_calls[node.call] += 1
        except KeyError:
            duplicate_calls[node.call] = 1
    duplicate_calls = [k for k in duplicate_calls if duplicate_calls[k] > 1]

    # Go through duplicates
    for duplicate_call in duplicate_calls:

        # Do this until we no longer are finding and merging neighbors
        while True:

            # This will control whether we break out of loop
            found_neighbor = False

            # Look for nodes whose call is duplicate_call
            these_dups = tree.search_nodes(**{"call":duplicate_call})

            # If we found such nodes
            if len(these_dups) > 1:

                # Go through the duplicates
                for dup in these_dups:

                    # Get clade containing dup and it's closest neighbor
                    anc = dup.get_ancestors()[0]
                    tips = anc.get_descendants()

                    # No single neighbor possible
                    if len(tips) != 2:
                        break

                    # If this nearest neighbor has the same call as the duplicate
                    if tips[0].call == tips[1].call:

                        # Figure out which is the dup and which is the neighbor
                        if tips[0].name == dup.name:
                            neighbor = tips[1]
                        else:
                            neighbor = tips[0]

                        # Record that are merging the neighbor into the duplicate node,
                        # node features (which will include name, call, etc.)
                        neigh_dict = {}
                        for f in neighbor.features:
                            try:
                                neigh_dict[f] = copy.deepcopy(neighbor.__dict__[f])
                            except KeyError:
                                pass

                        dup.merge_stack.append(neigh_dict)

                        # Delete the neighbor and break out of the loop
                        neighbor.delete()
                        found_neighbor = True
                        break

            # If we did not find a neighbor duplicate, break out of search
            if not found_neighbor:
                break

    return tree

def expand_sisters(some_tree,work_on_copy=True):
    """
    Reverse operation of merge_sisters.

    some_tree: ete3 tree with duplicate sisters merged using merge_sisters
    work_on_copy: work on a copy of the tree rather than the input tree

    returns tree with sister duplicates expanded
    """

    # Copy tree -- do not operate on input tree directly
    if work_on_copy:
        tree = some_tree.copy(method="deepcopy")
    else:
        tree = some_tree

    while True:

        found_expansion = False

        for node in tree.get_leaves():

            if len(node.merge_stack) > 0:

                # Pull nodes from merge_stack, reversing them so we access
                # as last-in, first-out.
                to_merge = node.merge_stack[::-1]
                node.merge_stack = []

                # Merge from last to first
                for m in to_merge:

                    # Add two children to node. One will be "node" duplicated,
                    # the other will be the merged sister.
                    new_node_1 = node.add_child()
                    new_node_2 = node.add_child()

                    # Copy the node features to new_node_1
                    new_node_1 = copy_node_features(node,new_node_1)

                    # Copy the features out of the merge dict into new_node_2.
                    for feature_key in m:
                        new_node_2.add_feature(feature_key,m[feature_key])

                    # new_node_1 is now the source for new additions from
                    # merge_stack.
                    node = new_node_1

                found_expansion = True
                break

        if not found_expansion:
            break

    return tree


def leaves_bifurcating_on_ancestor(some_tree,node_1,node_2):

    anc_tree = some_tree.get_common_ancestor([node_1,node_2])

    # Pre-order traverses root, left, right
    root_node = None
    left_anc_node = None
    for node in anc_tree.traverse("preorder"):

        # First iteration gets root node
        if root_node is None:
            root_node = node
            continue

        # Second iteration gets ancestor of all left. Get left leaves
        if left_anc_node is None:
            left_anc_node = node
            left_leaves = left_anc_node.get_leaves()
            continue

        # first node without left_anc_node as a descendant is the
        # right ancestor.
        if left_anc_node not in node.get_ancestors():
            right_anc_node = node
            right_leaves = right_anc_node.get_leaves()
            break

    out_dict = {"anc_node":root_node,
                "left_node":left_anc_node,
                "right_node":right_anc_node,
                "left_leaves":left_leaves,
                "right_leaves":right_leaves}

    return out_dict
