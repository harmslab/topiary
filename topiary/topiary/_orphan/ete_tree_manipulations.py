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

def build_species_corrected_gene_tree(df,species_tree,gene_tree_string):
    """
    Hacked attempt to build a species-corrected gene tree. Use at your peril.

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

def _copy_root(unrooted_newick,
               rooted_newick,
               output_newick,
               unrooted_tree_fmt=0,
               rooted_tree_fmt=0):
    """
    Root the tree in an unrooted newick file using the root from a rooted
    newick file with the same taxa. Write to an output file.

    unrooted_newick: newick file containing an unrooted tree
    rooted_newick: newick file containing a rooted tree with the same taxa
                   as the unrooted tree
    output_newick: output file to write results
    unrooted_tree_fmt: what to preserve from unrooted tree. integer.
                       interpretation is done by ETE3 (table below
                       current as of v. 3.1.1).
    rooted_tree_fmt:   what to preserve from rooted tree. integer.
                       interpretation is done by ETE3 (table below
                       current as of v. 3.1.1).

     |        ======  ==============================================
     |        FORMAT  DESCRIPTION
     |        ======  ==============================================
     |        0        flexible with support values
     |        1        flexible with internal node names
     |        2        all branches + leaf names + internal supports
     |        3        all branches + all names
     |        4        leaf branches + leaf names
     |        5        internal and leaf branches + leaf names
     |        6        internal branches + leaf names
     |        7        leaf branches + all names
     |        8        all names
     |        9        leaf names
     |        100      topology only
     |        ======  ==============================================

    """

    # Load trees
    rooted_tree = Tree(rooted_newick,format=rooted_tree_fmt)
    unrooted_tree = Tree(unrooted_newick,format=unrooted_tree_fmt)

    # Make sure they have the same taxa
    rooted_leaves = set(rooted_tree.get_leaf_names())
    unrooted_leaves = set(unrooted_tree.get_leaf_names())
    if rooted_leaves != unrooted_leaves:
        err = "both trees must have the exact same leaves\n"
        raise ValueError(err)

    left_leaves = []
    right_leaves = []

    root_node = None
    left_anc_node = None

    # Pre-order traverses root, left, right
    for node in rooted_tree.traverse("preorder"):

        # First iteration gets root node
        if root_node is None:
            root_node = node
            continue

        # Second iteration gets ancestor of all left. Get leaves.
        if left_anc_node is None:
            left_anc_node = node
            left_leaves = [l.name for l in left_anc_node.get_leaves()]
            continue

        # First node without left_anc_node as a descendant is the
        # right ancestor. Get leaves.
        if left_anc_node not in node.get_ancestors():
            right_anc_node = node
            right_leaves = right_anc_node.get_leaves()
            break

    # If we have single outgroups on the left or right, root on that
    if len(left_leaves) == 1:
        unrooted_tree.set_outgroup(left_leaves[0])
    elif len(right_leaves) == 1:
        unrooted_tree.set_outgroup(right_leaves[0])

    # Otherwise, try to root on last common ancestor of left or right. This
    # may throw error, depending on tree topology, so we try left first,
    # and then try right if that does not work.
    else:

        root_successful = False
        try:
            root_anc = unrooted_tree.get_common_ancestor(*left_leaves)
            unrooted_tree.set_outgroup(root_anc)
            root_successful = True
        except ete3.coretype.tree.TreeError:
            try:
                root_anc = unrooted_tree.get_common_ancestor(*right_leaves)
                unrooted_tree.set_outgroup(root_anc)
                root_successful = True
            except ete3.coretype.tree.TreeError:
                pass

        if not root_successful:
            unrooted_tree.set_outgroup(root_anc.get_children()[1])

    # Write out newly rooted tree
    unrooted_tree.write(outfile=output_newick)
