"""
Load a tree into an ete3 tree data structure.
"""

from topiary._private.check import check_bool

import ete3
from ete3 import Tree
import dendropy as dp

import glob
import os
import re

def read_tree(tree,fmt=None):
    """
    Load a tree into an ete3 tree data structure.

    Parameters
    ----------
    tree : ete3.Tree or dendropy.Tree or str
        some sort of tree. can be an ete3.Tree (returns self), a dendropy Tree
        (converts to newick and drops root), a newick file or a newick string.
    fmt : int or None
        format for reading tree from newick. 0-9 or 100. (See Notes for what
        these mean). If fmt is None, try to parse without a format descriptor,
        then these formats in numerical order.

    Returns
    -------
    tree : ete3.Tree
        an ete3 tree object.

    Notes
    -----
    `fmt` number is read directly by ete3. See their documentation for how these
    are read (http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees).
    As of ETE3.1.1, these numbers mean:

    + 0: flexible with support values
    + 1: flexible with internal node names
    + 2: all branches + leaf names + internal supports
    + 3: all branches + all names
    + 4: leaf branches + leaf names
    + 5: internal and leaf branches + leaf names
    + 6: internal branches + leaf names
    + 7: leaf branches + all names
    + 8: all names
    + 9: leaf names
    + 100: topology only

    """

    # Already an ete3 tree.
    if issubclass(type(tree),ete3.TreeNode):
        return tree

    # Convert dendropy tree into newick (drop root)
    if issubclass(type(tree),dp.Tree):
        tree = tree.as_string(schema="newick",suppress_rooting=True)

    # If we get here, we need to convert. If fmt is not specified, try to parse
    # without a format string.
    if fmt is None:

        try:
            t = Tree(tree)
        except ete3.parser.newick.NewickError:

            # Try all possible formats now, in succession
            w = "\n\nCould not parse tree without format string. Going to try different\n"
            w += "formats. Please check output carefully.\n\n"
            print(w)

            formats = list(range(10))
            formats.append(100)

            t = None
            for f in formats:
                try:
                    t = Tree(tree,format=f)
                    w = f"\n\nSuccessfully parsed tree with format style {f}.\n"
                    w += "Please see ete3 documentation for details:\n\n"
                    w += "http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees\n\n"
                    print(w)
                    break

                except ete3.parser.newick.NewickError:
                    continue

            if t is None:
                err = "\n\nCould not parse tree!\n\n"
                raise ValueError(err)

    else:
        # Try a conversion with the specified format
        t = Tree(tree,format=fmt)

    return t


def _map_tree_to_tree(T1,T2):
    """
    Map nodes between one tree and another based on shared descendants.

    Parameters
    ----------
    T1 : ete3.Tree or toytree.tree
        one tree to compare
    T2 : ete3.Tree or toytree.tree
        second tree to compare

    Returns
    -------
    shared_nodes : list
        list of tuples containing shared nodes between the two trees
    T1_nodes : list
        list of nodes from T1 that are not in T2
    T2_nodes : list
        list of nodes from T2 that are not in T1
    """

    def _ete3_node_dict(T):
        """
        Create dictionary keying ete3 tree nodes to tuple of descendants.
        """

        node_dict = {}
        for node in T.traverse():
            leaves = node.get_leaf_names()
            leaves = [t for t in leaves]
            leaves.sort()
            leaves = tuple(leaves)
            node_dict[leaves] = node

        return node_dict

    def _toytree_node_dict(T):
        """
        Create dictionary keying toytree tree nodes to tuple of descendants.
        """

        node_dict = {}
        for k in T.idx_dict:
            leaves = T.idx_dict[k].get_leaf_names()
            leaves.sort()
            leaves = tuple(leaves)
            node_dict[leaves] = T.idx_dict[k]

        return node_dict

    # Construct dictionary keying node to tuple of descendants for T1
    if issubclass(type(T1),ete3.Tree):
        T1_node_dict = _ete3_node_dict(T1)
    else:
        T1_node_dict = _toytree_node_dict(T1)

    # Construct dictionary keying node to tuple of descendants for T2
    if issubclass(type(T2),ete3.Tree):
        T2_node_dict = _ete3_node_dict(T2)
    else:
        T2_node_dict = _toytree_node_dict(T2)

    T1_nodes = []
    shared_nodes = []
    T2_keys = list(T2_node_dict.keys())
    for key in T1_node_dict:

        try:
            # Shared, not in T2 alone
            shared_nodes.append((T1_node_dict[key],T2_node_dict[key]))
            T2_keys.remove(key)

        except KeyError:
            # Only in T1
            T1_nodes.append(T1_node_dict[key])

    # If key is left in T2_keys, the node is only in T2
    T2_nodes = []
    for key in T2_keys:
        T2_nodes.append(T2_node_dict[key])


    return shared_nodes, T1_nodes, T2_nodes

def load_trees(directory=None,
               prefix=None,
               T_clean=None,
               T_support=None,
               T_anc_label=None,
               T_anc_pp=None,
               T_event=None):
    """
    Generate an ete3 tree with features 'event', 'anc_pp', 'anc_label',
    and 'bs_support' on internal nodes. This information is read from the input
    ete3 trees or the specified topiary output directory. The tree is rooted
    using T_event. If this tree is not specified, the midpoint root is used.
    Trees are read from the directory first, followed by any ete3 trees
    specified as arguments. (This allows the user to override trees from the
    directory if desired). If no trees are passed in, returns None.

    Warning: this will modify input ete3 trees as it works on the trees
    rather than copies.

    Parameters
    ----------
    directory : str
        output directory from a topiary calculation that has .newick files
        in it. Function will load all trees in that directory.
    prefix : str, optional
        what type of trees to plot from the directory. should be "reconciled"
        or "gene". If None, looks for reconciled trees. If it finds any, these
        prefix = "reconciled"
    T_clean : ete3.Tree, optional
        clean tree (leaf labels and branch lengths, nothing else). Stored as
        {}-tree.newick in output directories.
    T_support : ete3.Tree, optional
        support tree (leaf labels, branch lengths, supports). Stored as
        {}-tree_supports.newick in output directories.
    T_anc_label : ete3.Tree, optional
        ancestor label tree (leaf labels, branch lengths, internal names)
        Stored as {}-tree_anc-label.newick.
    T_anc_pp : ete3.Tree, optional
        ancestor posterior probability tree (leaf labels, branch lengths,
        posterior probabilities as supports) Stored as {}-tree_anc-pp.newick.
    T_event : ete3.Tree, optional
        tree with reconciliation events as internal labels (leaf labels,
        branch lengths, event labels). Stored as reconciled-tree_events.newick

    Returns
    -------
    merged_tree : ete3.Tree or None
        rooted tree with features on internal nodes. Return None if no trees
        are passed in.
    """

    # Load trees from the directory
    if directory is not None:

        tree_files = glob.glob(os.path.join(directory,"*.newick"))
        to_path = dict([(os.path.split(t)[-1],t) for t in tree_files])

        if prefix is None:

            num_rec = len([t for t in tree_files if t.startswith("reconciled")])
            if num_rec > 0:
                prefix = "reconciled"
            else:
                prefix = "gene"

        if T_clean is None:
            try:
                T_clean = ete3.Tree(to_path[f"{prefix}-tree.newick"],format=0)
            except KeyError:
                pass

        if T_support is None:
            try:
                T_support = ete3.Tree(to_path[f"{prefix}-tree_supports.newick"],format=0)
            except KeyError:
                pass

        if T_event is None:
            try:
                T_event = ete3.Tree(to_path[f"{prefix}-tree_events.newick"],format=1)
            except KeyError:
                pass

        if T_anc_label is None:
            try:
                T_anc_label = ete3.Tree(to_path[f"{prefix}-tree_anc-label.newick"],format=1)
            except KeyError:
                pass

        if T_anc_pp is None:
            try:
                T_anc_pp = ete3.Tree(to_path[f"{prefix}-tree_anc-pp.newick"],format=0)
            except KeyError:
                pass

    # This is the order of priority for getting branch lengths from the tree.
    # The output tree will be copied from the first non-None tree in this
    # list.
    T_list = [T_event,T_anc_pp,T_anc_label,T_support,T_clean]

    # Make sure trees were actually passed in. If none were, return None.
    T_list = [T for T in T_list if T is not None]
    if len(T_list) == 0:
        return None

    # Make sure all trees have the same descendants
    ref_leaves = set(list(T_list[0].get_leaf_names()))
    for T in T_list[1:]:
        test_leaves = set(list(T.get_leaf_names()))
        if len(set(ref_leaves) - set(test_leaves)) != 0:
            err = "All trees must have the same leaves.\n"
            raise ValueError(err)

    # If we have an event tree, root all trees on that rooted tree
    if prefix == "reconciled":

        if T_event is not None:
            root_tree = T_event
        else:
            root_tree = T_clean

        # Get left and right descendants of the root node
        root_on = []
        for n in root_tree.get_tree_root().iter_descendants():
            leaves = n.get_leaf_names()
            leaves.sort()
            root_on.append(tuple(leaves))
            if len(root_on) == 2:
                break

    # If not a reconciled tree, do midpoint rooting using first tree in list.
    else:
        root_on = [T_list[0].get_midpoint_outgroup().get_leaf_names()]
        root_on.append(list(set(T_list[0].get_leaf_names()) - set(root_on[0])))

    # For each tree...
    for T in T_list:

        # Do not root the event tree as it will wipe out label
        if T is T_event:
            continue

        # Get MRCA for descendants of left and right from event root
        T.unroot()
        left = T.get_common_ancestor(root_on[0])
        right = T.get_common_ancestor(root_on[1])

        # If there is only one descendant on one lineage, left and right
        # will be the exact same node. Set the single descendant, rather
        # than the MRCA as the outgroup node.
        if left is right:
            if len(root_on[0]) == 1:
                left = root_on[0][0]
            elif len(root_on[1]) == 1:
                right = root_on[1][0]
            else:
                w = "Should not have gotten here. Weird rooting state possible."
                print(w)

        # Try setting left, then right outgroup.
        try:
            T.set_outgroup(left)
        except ete3.coretype.tree.TreeError:
            T.set_outgroup(right)

    # Make new tree from first tree in list. This will be our output tree.
    out_tree = T_list[0].copy()
    for n in out_tree.traverse():

        # Create empty features
        if not n.is_leaf():
            n.add_feature("event",None)
            n.add_feature("anc_pp",None)
            n.add_feature("anc_label",None)
            n.add_feature("bs_support",None)


    # features_to_load maps trees with information to copy (keys) to what
    # feature we should extract from that tree. Values are (feature_in_tree,
    # name_of_feature_in_out_tree,whether_to_allow_root_value).
    features_to_load = {T_event:("name","event",True),
                        T_anc_pp:("support","anc_pp",False),
                        T_anc_label:("name","anc_label",False),
                        T_support:("support","bs_support",False)}

    # Only copy from trees that are not None.
    trees = [k for k in features_to_load.keys() if k is not None]

    stash_values = {}
    for T in trees:

        in_feature = features_to_load[T][0]
        out_feature = features_to_load[T][1]
        root_allowed = features_to_load[T][2]

        # Since all have the same root and descendants, tree nodes should be
        # identical between trees and uniquely identified by their descendants
        shared, T_alone, out_alone = _map_tree_to_tree(T,out_tree)
        if len(T_alone) > 0 or len(out_alone) > 0:
            err = "Cannot merge trees with different topologies.\n"
            raise ValueError(err)

        # Map data from features on input tree to the features on the output
        # tree
        for s in shared:

            # Node to copy from to
            in_node = s[0]
            out_node = s[1]

            # If internal
            if not out_node.is_leaf():

                # Get value from in
                try:
                    value = in_node.__dict__[in_feature]
                except KeyError:
                    value = in_node.__dict__[f"_{in_feature}"]

                # small hack --> anc to a
                if T is T_anc_label:
                    value = re.sub("anc","a",value)

                # Add value to out
                out_node.add_feature(out_feature,value)

            # If root node, pull out anc_ if there.
            if out_node.is_root():
                if not root_allowed:
                    stash_values[out_feature] = value
                    out_node.add_feature(out_feature,None)


    # Copy ancestor to correct node because displayed by rooting
    if len(stash_values) > 0 and "anc_label" in stash_values:

        if stash_values["anc_label"] != "":

            for n in out_tree.traverse():
                if n.is_leaf():
                    continue

                if n.anc_label == "":
                    n.add_feature("anc_label",stash_values["anc_label"])
                    n.add_feature("anc_pp",stash_values["anc_pp"])

    return out_tree

def write_trees(T,
                name_dict=None,
                out_file=None,
                overwrite=False,
                anc_pp=True,
                anc_label=True,
                bs_support=True,
                event=True):
    """
    Write out an ete3.Tree as a newick format. This function looks for features
    set by :code:`load_trees` and then writes an individual tree out with each
    feature. The features are :code:`anc_pp`, :code:`anc_label`, :code:`bs_support`,
    and :event:`event`. This will write out trees for any of these features 
    present; not all features need to be in place for this function to work. 

    Parameters
    ----------
    T : ete3.TreeNode
        ete3 tree with information loaded into appropriate features. This is the
        tree returned by :code:`load_trees`. 
    name_dict : dict
        name_dict : dict, optional
        dictionary mapping strings in node.name to more useful names. (Can be
        generated using :code:`topiary.draw.core.create_name_dict`). If not 
        specified, trees are written out with uid as tip names
    out_file : str, optional
        output file. If defined, write the newick string the file.
    overwrite : bool, default=False
        whether or not to overwrite an existing file
    anc_pp : bool, default=True
        whether or not to write a tree with anc_pp as support values
    anc_label : bool, default=True
        whether or not to write a tree with anc_label as internal node names
    bs_support : bool, default=True
        whether or not to write a tree with bs_support as support values
    event : bool, default=True
        whether or not to write a tree with events as internal node names

    Returns
    -------
    tree : str
        Newick string representation of the output tree(s)
    """
    
    # --------------------------------------------------------------------------
    # Parameter sanity checking

    if not issubclass(type(T),ete3.TreeNode):
        err = "\nT must be an ete3.Tree instance\n\n"
        raise ValueError(err)

    if name_dict is not None:
        if not issubclass(type(name_dict),dict):
            err = "\nname_dict must be a dictionary\n\n"
            raise ValueError(err)
    
    if out_file is not None:
        if not issubclass(type(out_file),str):
            err = "\nout_file must be a string pointing to a file to write out\n\n"
            raise ValueError(err)
        
        overwrite = check_bool(overwrite,"overwrite")

        if os.path.exists(out_file):
            if os.path.isfile(out_file):
                if overwrite:
                    os.remove(out_file)
                else:
                    err = f"\nout_file '{out_file}' exists. Either delete or set overwrite to True\n\n"
                    raise FileExistsError(err)
            else:
                err = f"\nout_file '{out_file}' exists but is a directory. Cannot write output.\n\n"
                raise FileExistsError(err)

    anc_pp = check_bool(anc_pp,"anc_pp")
    anc_label = check_bool(anc_label,"anc_label")
    bs_support = check_bool(bs_support,"bs_support")
    event = check_bool(event,"event")

    # --------------------------------------------------------------------------
    # Set up output

    out_trees = []

    # Work on a copy
    T = T.copy()

    # If name dict is specified
    if name_dict is not None:
        for n in T.traverse():
            if n.is_leaf():
                n.name = name_dict[n.name]

    # --------------------------------------------------------------------------
    # Create bs_supports, events, anc_label, and anc_pp

    # bs_supports
    if bs_support:
        write_bs_supports = False
        for n in T.traverse():
            if not n.is_leaf():
                if n.bs_support is not None:
                    n.support = n.bs_support
                    write_bs_supports = True

        # format=2 --> all branches + leaf names + internal supports
        if write_bs_supports:
            out_trees.append(T.write(format=2))

    # event
    if event:
        write_events = False
        for n in T.traverse():
            if not n.is_leaf():
                if n.event is not None:
                    n.name = n.event
                    write_events = True

        # format=3 --> all branches + all names
        if write_events:
            out_trees.append(T.write(format=3))

    # anc_label
    if anc_label:
        write_anc_label = False
        for n in T.traverse():
            if not n.is_leaf():
                if n.anc_label is not None:
                    n.name = n.anc_label
                    write_anc_label = True

        # format=3 --> all branches + all names
        if write_anc_label:
            out_trees.append(T.write(format=3))

    # anc_pp
    if anc_pp:
        write_anc_pp = False
        for n in T.traverse():
            if not n.is_leaf():
                if n.anc_pp is not None:
                    n.support = n.anc_pp
                    write_anc_pp = True

        # format=2 --> all branches + leaf names + internal supports
        if write_anc_pp:
            out_trees.append(T.write(format=2))

    # --------------------------------------------------------------------------
    # Finalize output

    out_trees = "\n".join(out_trees)

    if out_file is not None:
        f = open(out_file,'w')
        f.write(out_trees)
        f.close()

    return out_trees 
    