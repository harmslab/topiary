"""
Functions and classes for plotting toytree trees.
"""

import topiary
import ete3
import toytree
import toyplot

from toyplot import pdf, svg, png

import numpy as np
import re, os, copy

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

def _color_to_css(color):
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
        color = _color_to_css(color)
        cm = lambda prop: color
        span = []

    # color dictionary
    elif issubclass(type(color),dict):

        prop_keys = list(set(prop))
        for k in prop_keys:
            try:
                color[k] = _color_to_css(color[k])
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
            color = _color_to_css(color)
            cm = lambda prop: color
            span = []

        # color gradient
        elif len(color) == 2:

            # Make sure the property is a float
            try:
                prop = np.array(prop,dtype=float)
                if np.sum(np.isnan(prop)) > 0:
                    raise ValueError
            except (ValueError,TypeError):
                err = "\nProperty is not a float. Using a color\n"
                err += "gradient requires a float property value."
                raise ValueError(err)

            # Create palette if needed
            if palette is None:
                if color is not None:
                    color = list(color)
                    for i in range(len(color)):
                        color[i] = _color_to_css(color[i])

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



def load_trees(prev_run_dir,
               tree_base_path,
               tree_names,
               tree_fmt=None,
               outgroup=None):
    """
    Load trees and tree formatting information from a previous run directory.

    Parameters
    ----------
    prev_run_dir : str
        previous run directory (for example, output from an ML tree inference)
    tree_base_path : str
        base directory to look for trees within prev_run_dir. (for example, look
        in ml_tree/output for output trees for ML tree inference)
    tree_names : list
        tree files to look for in tree_base_path
    tree_fmt : list, default=1
        ete3 Tree format codes to use for reading trees. Should be ints
        indiciating which tree format should be applied to read each file in
        tree_names. Should be the same length as tree_names.
    df : pandas.DataFramem, optional
        topiary dataframe to use to interpret tip_columns. overrides the
        dataframe that is pulled from prev_run_dir.

    Return
    ------
    tree : ete3.Tree or list
        tree instances for each tree requested in tree_names.
    """

    # If the user did not specify tree formats, make them all 1
    if tree_fmt is None:
        tree_fmt = [1 for _ in tree_names]

    # Check tree format and tree names compatibility
    if len(tree_names) != len(tree_fmt):
        err = "\ntree_fmt length does not match tree_names length\n\n"
        raise ValueError(err)

    # Go over all tree names
    T_list = []
    for i, t in enumerate(tree_names):

        tree_path = os.path.join(prev_run_dir,tree_base_path,t)
        if not os.path.exists(tree_path):
            err = f"\nCould not find tree file '{tree_path}'\n\n"
            raise FileNotFoundError(err)

        # read tree
        T = topiary.io.read_tree(tree_path,fmt=tree_fmt[i])

        # If outgroup specified
        if outgroup is not None:

            # Root tree. We have to try with both the left and right
            # descendants of the previous root.
            try:
                common_anc = T.get_common_ancestor(outgroup[0])
                T.set_outgroup(common_anc)
            except ete3.coretype.tree.TreeError:
                common_anc = T.get_common_ancestor(outgroup[1])
                T.set_outgroup(common_anc)

        T_list.append(T)

    if len(T_list) == 1:
        return T_list[0]

    return T_list

def create_name_dict(df,tip_columns=None,separator="|"):
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
    for name in name_to_uid:
        if len(name_to_uid[name]) > 1:
            for uid in name_to_uid[name]:
                uid_to_name[uid] = separator.join([uid_to_name[uid],uid])

    return uid_to_name


def final_render(pt,output_file=None,default_file="tree.pdf"):
    """
    Render tree on screen and/or file. Decides what file type to write to based
    on extension (svg,png,pdf). This function is a a kludge, but exists so the
    user-facing functions (like ml_tree) can have output_file=None as their
    default without making those functions aware of whether this is running in
    a notebook or not.

    Parameters
    ----------
    pt : prettyTree
        object to render
    output_file : str, optional
        write tree to specified output_file. If None and in a notebook,
        just draw on the screen. If None and *not* in a notebook,
        write only to default_file
    default_file : str, default="tree.pdf"
        write to this file. if None, do not write.

    Returns
    -------
    canvas : toytree.canvas or None
        toytree.canvas (if in a notebook) or None (if not in a notebook)
    """

    # If output_file is not specified...
    if output_file is None:

        # If not in a notebook, create a file
        if not topiary._in_notebook:
            output_file = default_file

    # If output_file specified, write out
    if output_file is not None:
        key = output_file[-3:].lower()
        render_dict = {"svg":toyplot.svg.render,
                       "pdf":toyplot.pdf.render,
                       "png":toyplot.png.render}
        try:
            render_fcn = render_dict[key]
        except KeyError:
            print(f"Could not identify render type for file {output_file}.")
            print("Using pdf.",flush=True)
            render_fcn = render_dict["pdf"]

        render_fcn(pt.canvas,output_file)

    # If this is in a notebook, render it and return so it appears
    if topiary._in_notebook:
        return pt.canvas

    return None
