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
        :code:`["species","recip_paralog"], then :code:`["species","nickname"]`,
        then :code:`["species",f"{name[:10]}..."]`.
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
        if "recip_paralog" in df.columns:
            tip_columns = ["species","recip_paralog"]
        elif "nickname" in df.columns:
            tip_columns = ["species","nickname"]
        else:
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
    span : list
        list of values corresponding to span. For continuous map, will be 
        [min,max]; for discrete map will be values for each color
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
