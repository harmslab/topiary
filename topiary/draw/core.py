"""
Functions and classes for plotting toytree trees.
"""

import topiary
import ete3
import toytree
import toyplot

import numpy as np
import re
import os
import copy
import glob

def _protect_name(name):
    """
    Protect a name going into a newick string.
    """

    new_name = re.sub(" ","%20",name)
    return f"'{new_name}'"

def _deprotect_name(name):
    """
    Invert _protect_name for bringing name out of newick string.
    """

    return re.sub("%20"," ",name).strip("'")

def color_to_css(color):
    """
    Take a color (as hex string, rgb tuple, rgba tuple, named color, or
    toyplot css) and return toyplot css.
    """

    if issubclass(type(color),str):
        color = toyplot.color.broadcast(color,1)[0]

    elif hasattr(color,"__iter__"):
        if len(color) == 3:
            color = toyplot.color.rgb(*color)
        elif len(color) == 4:
            color = toyplot.color.rgba(*color)
        else:
            err = f"color '{color}' not recognized. Should be a named color,\n"
            err += "RGB tuple, RGBA tuple, or css.\n\n"
            raise ValueError(err)

    else:
        err = f"color '{color}' not recognized. Should be a named color,\n"
        err += "RGB tuple, RGBA tuple, or css.\n\n"
        raise ValueError(err)

    return toyplot.color.to_css(color)


def get_round_to(value,total_requested=2):
    """
    Figure out a natural place to round a float for creating useful but pretty
    strings.

    Parameters
    ----------
    value : float
        value to round
    total_requested : int, default=3
        target number of digits for rounding. If value is 1234.23 and
        total_requested is 3, this would round to 1234. If total_requested
        is 5, this would give 1234.23.

    Returns
    -------
    round_to : int
        what number to round value to prior to string conversion
    """

    # Convert value to string, dropping sign
    value_str = f"{np.abs(value)}"

    # Deal with exponents (1e+50, for example)
    if re.search("e",value_str):

        # Split on "e"
        exp_split = value_str.split("e")
        base = exp_split[0]
        exponent = int(exp_split[1])

        # Process base
        base_split = base.split(".")
        base_whole = list(base_split[0])
        if len(base_split) == 2:
            base_decimal = list(base_split[1])
        else:
            base_decimal = []

        # Positive exponent
        if exponent > 0:
            decimal_out = ["0" for _ in range(int(exponent))]
            grab = min([len(decimal_out),len(base_decimal)])
            decimal_out[:grab] = base_decimal[:grab]
            value_str = "".join(base_whole) + "".join(decimal_out)

        # Negative exponent
        else:
            decimal_out = ["0" for _ in range(-exponent)]
            grab = min([len(decimal_out),len(base_decimal)])
            for i in range(grab):
                decimal_out[-(i+1)] = base_decimal[i]
            decimal_out[-(grab + len(base_whole)):-grab] = base_whole[:]

            value_str = "0." + "".join(decimal_out)

    # Split value on "."
    value_split = f"{value_str}".split(".")

    # If there is no decimal, round at 0
    if len(value_split) == 1:
        round_at = 0
    else:

        # Consider whole and decimal portions of the float
        whole = value_split[0]
        decimal = value_split[1]

        non_zero_num_whole = len(whole.lstrip("0"))

        round_at = 0
        num_taken = non_zero_num_whole
        for counter, d in enumerate(decimal):

            # If we've seen at least one non-zero digit and we've reached
            # requested total digits, break
            if num_taken >= total_requested:
                if non_zero_num_whole > 0 or round_at > 0:
                    break

            num_taken += 1
            if d != "0":
                round_at = (counter + 1)

    return round_at


def create_name_dict(df,tip_columns=None,separator="|",disambiguate=True):
    """
    Create a dictionary mapping between uid and pretty names extracted from
    tip_columns in the dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe
    tip_columns : list, optional
        columns in dataframe to use to create human-readable names for tips on
        tree. For example, :code:`tip_columns = ["species","nickname"]` would
        give names like :code:`"Homo sapiens|LY96"`. If None, try
        :code:`["species","nickname"]`. If nickname is not in dataframe,
        fall back to :code:`["species",f"{name[:10]}..."]`.
    separator : str
        separator to use between pretty names. Cannot be "#,;:'\")(" as these
        are used in newick format.
    disambiguate : bool, default=True
        if two tip labels will be the same (for example, two labels will be
        "Homo sapiens|LY96"), append the uid to those labels so they can be
        uniquely identified.

    Return
    ------
    name_dict : dict
        dict mapping between uid and pretty name
    """

    # Check sanity of tip_name_separator
    if separator in topiary._private.reserved_characters:
        err = f"\nseparator cannot be {separator}. Disallowed values are:\n\n"
        for r in topiary._private.reserved_characters:
            err += f"    {r}\n"
        err += "\n\n"
        raise ValueError(err)

    local_df = df.loc[df.keep,:]
    if tip_columns is None:

        try:
            local_df.nickname
            tip_columns = ["species","nickname"]
        except AttributeError:

            local_df = df.copy()
            trunc_name = []
            for i in range(len(local_df)):
                if len(local_df["name"].iloc[i]) > 10:
                    trunc_name.append(f"{local_df['name'].iloc[i][:10]}...")
                else:
                    trunc_name.append(f"{local_df['name'].iloc[i]}")

            local_df["trunc_name"] = trunc_name
            tip_columns = ["species","trunc_name"]

    # Make sure the tip columns are sane
    for c in tip_columns:
        try:
            local_df[c]
        except KeyError:
            err = f"\ncolumn '{c}' not in dataframe. Please check the\n"
            err += "tip_columns argument.\n\n"
            raise ValueError(err)

    # Construct uid_to_name, which will map between the uid on the leaves to a
    # prettier/more useful name. Store both uid_to_name and name_to_uid. The
    # second dictionary allows us to look for duplicated names.
    uid_to_name = {}
    name_to_uid = {}
    for i in range(len(local_df)):

        uid = local_df.uid.iloc[i]

        name = []
        for column in tip_columns:
            name.append(f"{local_df[column].iloc[i]}")
        name = separator.join(name)

        uid_to_name[uid] = name

        try:
            name_to_uid[name].append(uid)
        except KeyError:
            name_to_uid[name] = [uid]

    # Look for duplicated names. Append uid to duplicated names to make unique.
    if disambiguate:
        for name in name_to_uid:
            if len(name_to_uid[name]) > 1:
                for uid in name_to_uid[name]:
                    uid_to_name[uid] = separator.join([uid_to_name[uid],uid])

    return uid_to_name


def map_tree_to_tree(T1,T2):
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


def ete3_to_toytree(T):
    """
    Generate a toytree.tree instance from an ete3.Tree instance. Copies over
    branch lengths, labels, and features.

    Parameters
    ----------
    T : ete3.Tree
        ete3.Tree to be converted

    Returns
    -------
    tT : toytree.tree
        converted toytree.tree
    """

    all_leaf_names = list(T.get_leaf_names())
    if len(set(all_leaf_names)) != len(all_leaf_names):
        err = "\nete3_to_toytree requires each leaf have a unique name.\n\n"
        raise ValueError(err)

    # Write out a tree with names protected (quotes, spaces removed, etc.)
    to_write_T = copy.deepcopy(T)
    for node in to_write_T.traverse():
        node.name = _protect_name(node.name)

    tT = toytree.tree(to_write_T.write())

    # Clean up names that toytree read in
    for k in tT.idx_dict:
        name = _deprotect_name(tT.idx_dict[k].name)
        tT.idx_dict[k].add_feature("name",name)

    # Map nodes unambiguously by descendant leaf names
    # toytree
    tT_node_dict = {}
    for k in tT.idx_dict:
        tT_leaves = tT.idx_dict[k].get_leaf_names()
        tT_leaves.sort()
        tT_leaves = tuple(tT_leaves)
        tT_node_dict[tT_leaves] = tT.idx_dict[k]

    # ete3.Tree
    T_node_dict = {}
    for node in T.traverse():
        T_leaves = node.get_leaf_names()
        T_leaves = [_deprotect_name(t) for t in T_leaves]
        T_leaves.sort()
        T_leaves = tuple(T_leaves)
        T_node_dict[T_leaves] = node

    # We can now map between toytree and ete3.Trees based on their shared keys.
    for node in tT_node_dict:

        # Get equivalent toytree and ete3 nodes
        tT_node = tT_node_dict[node]
        T_node = T_node_dict[node]

        # Copy all features from ete3 node to toytree node
        for f in T_node.features:

            # Manually added features and name will have format __dict__[key].
            # Reserved features (dist, support, minimally) have format
            # __dict__[_key].
            try:
                value = T_node.__dict__[f]
            except KeyError:
                value = T_node.__dict__[f"_{f}"]

            # Add as feature to toytree Tree
            tT_node.add_feature(f,value)

    return tT

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

        if T_event is None:
            err = "\nTree prefix 'reconciled' cannot be used without an event\n"
            err += "tree in the directory.\n\n"
            raise ValueError(err)

        # Get left and right descendants of the root node
        root_on = []
        for n in T_event.get_tree_root().iter_descendants():
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
        shared, T_alone, out_alone = map_tree_to_tree(T,out_tree)
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

def construct_colormap(color,prop,prop_span=None,palette=None):
    """
    Construct a toyplot color map. If value_min == value_max, the cmap will
    always return color_min.

    Parameters
    ----------
    color : str or tuple or dict
        set node color. If a single value, color all nodes that color. If
        list-like and length 2, treat as colors for minimum and maximum of a
        color gradient.  If dict, map property keys to color values. Colors
        can be RGBA tuples, named colors, or hexadecimal strings. See the
        toyplot documentation for allowed values.
    prop : list-like
        values of properties over which to make the color map
    prop_span : tuple, optional
        set min/max values for property color/size calculation. First element
        is min, second is max. Only applicable if color and size are set
        to be gradients.
    palette : toyplot.color.Palette, optional
        custom toyplot palette. if specified, use this palette rather than
        color to define color scheme for a color gradient. requires the
        property value be a float.

    Returns
    -------
    colormap : toyplot.LinearMap
        toyplot linear colormap
    """

    # This is a hack that makes sure the code below goes to the color
    # gradient calc, regardless of what the user passes in for color
    if palette is not None:
        color = (None,None)

    # color string/hexadecimal
    if issubclass(type(color),str):
        color = color_to_css(color)
        cm = lambda prop: color
        span = []

    # color dictionary
    elif issubclass(type(color),dict):

        prop_keys = list(set(prop))
        for k in prop_keys:
            try:
                color[k] = color_to_css(color[k])
            except KeyError:
                err = f"\nproperty value '{k}' not in color dictionary\n\n"
                raise ValueError(err)
        cm = lambda prop: color[prop]

        # Get keys that are seen in prop_keys. (Iterate in this way to keep
        # order the same as color.keys()).
        span = [k for k in color.keys() if k in prop_keys]


    # rgb, rgba, or gradient
    elif hasattr(color,"__iter__"):

        # rgb or rgba
        if len(color) in [3,4]:
            color = color_to_css(color)
            cm = lambda prop: color
            span = []

        # color gradient
        elif len(color) == 2:

            # Make sure the property is a float
            try:
                prop = np.array(prop,dtype=float)
            except (ValueError,TypeError):
                err = "\nProperty is not a float. Using a color\n"
                err += "gradient requires a float property value."
                raise ValueError(err)

            # Create palette if needed
            if palette is None:
                if color is not None:
                    color = list(color)
                    for i in range(len(color)):
                        color[i] = color_to_css(color[i])

                palette = toyplot.color.Palette(colors=color)

            # Check palette sanity
            if not issubclass(type(palette),toyplot.color.Palette):
                err = f"\nPallet {palette} should be an instance of toyplot.color.Palette\n\n"
                raise TypeError(err)

            if prop_span is None:
                min_value = np.min(prop)
                max_value = np.max(prop)
            else:
                min_value = prop_span[0]
                max_value = prop_span[1]

            # Construct colormap
            colormap = toyplot.color.LinearMap(palette=palette,
                                               domain_min=min_value,
                                               domain_max=max_value)
            cm = colormap.css
            if min_value == max_value:
                span = []
            else:
                span = [min_value,max_value]

        else:
            err = f"\ncolor '{color}' not recognized.\n\n"
            raise ValueError(err)

    else:
        err = f"\ncolor '{color}' not recognized.\n\n"
        raise ValueError(err)

    return cm, span

def construct_sizemap(size,prop,prop_span=None):
    """
    Construct a sizemap function that gives size given a property value.

    Parameters
    ----------
    size : float or tuple or dict, optional
        set node size. If a single value, make all nodes that size. If
        list-like and length 2, treat as sizes for minimum and maximum of a
        size gradient. If dict, map property keys to size values. Sizes must
        be float >= 0.
    prop : list-like
        values of properties over which to make the color map
    prop_span : tuple, optional
        set min/max values for property color/size calculation. First element
        is min, second is max. Only applicable if color and size are set
        to be gradients.

    Returns
    -------
    sm : function
        function that yields a size given the property
    """

    # single float
    try:
        size = float(size)
        if size < 0:
            err = "\nNode size must be >= 0\n\n"
            raise ValueError(err)
        sm = lambda prop: size
        span = []
        return sm, span

    except:
        pass

    # color dictionary
    if issubclass(type(size),dict):
        sm = lambda prop: size[prop]
        span = list(size.keys())
        return sm, span

    # size gradient
    elif hasattr(size,"__iter__"):

        # size gradient
        if len(size) == 2:

            # Parse size into two float
            size_min = topiary._private.check.check_float(size[0],
                                                          "size[0]",
                                                          minimum_allowed=0)
            size_max = topiary._private.check.check_float(size[1],
                                                          "size[1]",
                                                          minimum_allowed=0)

            # Make sure the property is a float
            try:
                prop = np.array(prop,dtype=float)
                if np.sum(np.isnan(prop)) > 0:
                    raise ValueError

                value_min = np.min(prop)
                value_max = np.max(prop)

            except (ValueError,TypeError):
                err = "\nProperty is not a float. Using a size\n"
                err += "gradient requires a float property value."
                raise ValueError(err)

            if prop_span is not None:
                value_min = prop_span[0]
                value_max = prop_span[1]

            # No slope in property or size -- just return minimum size
            if value_max == value_min or size_min == size_max:
                sm = lambda prop: size_min
                span = []
                return sm, span

            # Slope in property
            else:
                slope = (size_max - size_min)/(value_max - value_min)

                # Calculate intercept on size
                intercept = size_max - (slope*value_max)

                if slope == 0:
                    intercept = size_min

                sm = lambda prop: slope*prop + intercept
                span = [value_min,value_max]
                return sm, span

        else:
            err = f"\nsize '{size}' not recognized.\n\n"
            raise ValueError(err)

    else:
        err = f"\nsize '{size}' not recognized.\n\n"
        raise ValueError(err)


def parse_span_color(span_color,color):
    """
    Decide how to treat span_color (bs_color or pp_color), returning output that
    can feed into PrettyTree.draw_nodes().

    Parameters
    ----------
    span_color :
        putative bs_color or pp_color dictionary passed in by user
    color :
        putative color passed in by user

    Returns
    -------
    plot_nodes : bool
        whether or not we're plotting support nodes
    value_span : list or None
        top and bottom values for color map
    color_span : list or None
        top and bottom colors for color map. if "color" is specified, just
        return this.
    """

    # color takes precedence over span_color
    if color is not None:
        return False, None, color

    # if span_color is None,  return nothing. PrettyTree will do default color
    if span_color is None:
        return False, None, None

    bad_span_color = False
    if not issubclass(type(span_color),dict):
        bad_span_color = True
    else:

        keys = list(span_color.keys())
        try:
            min_value = float(keys[0])
            max_value = float(keys[1])
            if len(keys) != 2:
                raise IndexError
        except (KeyError,IndexError):
            bad_span_color = True

    if bad_span_color:
        err = "bs_color and pp_color must be dictionaries with two elements. The\n"
        err += "keys should be minimum and maximum support values for the color\n"
        err += "map. The values should be the bottom and top colors.\n"
        raise ValueError(err)

    span_span = (min_value,max_value)
    span_color = tuple(span_color.values())

    return True, span_span, span_color

def parse_position_string(position,x_offset,y_offset):
    """
    Get displacement in x and y for label placement. Displacement will be in
    tree coordinates.

    Parameters
    ----------
    position : str
        Allowed values are "top-left", "top", "top-right", "right",
        "bottom-right", "bottom", "bottom-left", and "left".
    x_offset : float
        how far to move in x
    t_offset : float
        how far to move in t

    Returns
    dx : float
        displacement in x
    dy : float
        displacement in y
    """

    axes = {"top":1,"bottom":1,"left":0,"right":0}
    moves = {"top":1,"bottom":-1,"left":-1,"right":1}
    entries = position.split("-")

    out = [0,0]
    for e in entries:
        try:
            axis = axes[e]
        except KeyError:
            err = f"position formatter '{e}' not recognized. Should be one of\n"
            err += "'top', 'bottom', 'left', 'right'.\n"
            raise ValueError(err)

        move = moves[e]
        out[axis] = move

    return out[0]*x_offset, out[1]*y_offset
