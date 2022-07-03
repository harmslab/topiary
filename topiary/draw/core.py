"""
Functions and classes for plotting ete3 trees.
"""

import topiary
import ete3

import numpy as np
import re, os

class ColorMap:
    """
    Class for mapping a value to a specific RGB value, returned as a hexadecimal
    value interpretable by ete3.
    """

    def __init__(self,span=(0,1),R=(1,0),G=(1,0),B=(1,0)):
        """
        Initialize class instance.

        Parameters
        ----------
            span: tuple holding minimum and maximum values for color map.
            R: tuple holding red channel value corresponding to minimum and maximum
               values of span. Values should be between 0 and 1.
            G: tuple holding green channel value corresponding to minimum and maximum
               values of span. Values should be between 0 and 1.
            B: tuple holding blue channel value corresponding to minimum and maximum
               values of span. Values should be between 0 and 1.
        """

        # Deal with span
        try:
            if len(span) != 2:
                raise ValueError

            self._min_value = float(span[0])
            self._max_value = float(span[1])

            if self._min_value >= self._max_value:
                raise ValueError

        except (IndexError,ValueError,TypeError,KeyError):
            err = "\nspan should be a length-2 tuple holding the minimum and\n"
            err += "maximum values corresponding to the low and high values \n"
            err += "of the ColorMap. For example, if bootstrap support values\n"
            err += "range from 0-100, span should be (0,100). The second value\n"
            err += "must be larger than the first value.\n"
            raise ValueError(err)

        # Deal with red channel
        try:
            if len(R) != 2:
                raise ValueError

            self._R_bottom = float(R[0])
            self._R_top = float(R[1])

            if self._R_bottom < 0 or self._R_bottom > 1:
                raise ValueError

            if self._R_top < 0 or self._R_top > 1:
                raise ValueError

        except (IndexError,ValueError,TypeError,KeyError):
            err = "\nR should be a length-2 tuple holding the minimum and\n"
            err += "maximum red channel values that correspond to the low and\n"
            err += "high values of this channel in the ColorMap. For example,\n"
            err += "if the channel should go from off to on over the map, this\n"
            err += "should be (0,1). These values must be between 0 and 1.\n"
            raise ValueError(err)

        # Determine slope and intercept for R channel
        self._R_slope = (self._R_top - self._R_bottom)
        self._R_intercept = self._R_bottom

        # Deal with green channel
        try:
            if len(G) != 2:
                raise ValueError

            self._G_bottom = float(G[0])
            self._G_top = float(G[1])

            if self._G_bottom < 0 or self._G_bottom > 1:
                raise ValueError

            if self._G_top < 0 or self._G_top > 1:
                raise ValueError

        except (IndexError,ValueError,TypeError,KeyError):
            err = "\nG should be a length-2 tuple holding the minimum and\n"
            err += "maximum green channel values that correspond to the low and\n"
            err += "high values of this channel in the ColorMap. For example,\n"
            err += "if the channel should go from off to on over the map, this\n"
            err += "should be (0,1). These values must be between 0 and 1.\n"
            raise ValueError(err)

        # Determine slope and intercept for R channel
        self._G_slope = (self._G_top - self._G_bottom)
        self._G_intercept = self._G_bottom

        # Deal with green channel
        try:
            if len(B) != 2:
                raise ValueError

            self._B_bottom = float(B[0])
            self._B_top = float(B[1])

            if self._B_bottom < 0 or self._B_bottom > 1:
                raise ValueError

            if self._B_top < 0 or self._B_top > 1:
                raise ValueError

        except (IndexError,ValueError,TypeError,KeyError):
            err = "\nGBshould be a length-2 tuple holding the minimum and\n"
            err += "maximum blue channel values that correspond to the low and\n"
            err += "high values of this channel in the ColorMap. For example,\n"
            err += "if the channel should go from off to on over the map, this\n"
            err += "should be (0,1). These values must be between 0 and 1.\n"
            raise ValueError(err)

        # Determine slope and intercept for R channel
        self._B_slope = (self._B_top - self._B_bottom)
        self._B_intercept = self._B_bottom


    def hex(self,value):
        """
        Return the hexadecimal representation of the RGB value corresponding to
        "value" given the span and R,G,B defined when the class was initialized.
        """

        try:
            value = float(value)
        except (ValueError,TypeError):
            err = "\nvalue must be a float\n\n"
            raise ValueError(err)

        # Force value to be within specified min/max range
        if value < self._min_value:
            value = self._min_value
        if value > self._max_value:
            value = self._max_value

        # Scale value to between 0 and 1.
        scaled = (value - self._min_value)/(self._max_value - self._min_value)

        # Determine RGB for this value
        R = int(round(255*(self._R_intercept + self._R_slope*scaled),0))
        G = int(round(255*(self._G_intercept + self._G_slope*scaled),0))
        B = int(round(255*(self._B_intercept + self._B_slope*scaled),0))

        # Return as a string HEX value
        return f"#{bytearray((R,G,B)).hex()}"

def load_trees(prev_run_dir,
               tree_base_path,
               tree_names,
               tree_fmt=None,
               outgroup=None):
    """
    Load trees and tree formatting information from a previous run directory.

    Parameters
    ----------
        prev_run_dir: previous run directory (for example, output from an ML tree
                      inference)
        tree_base_path: base directory to look for trees within prev_run_dir. (for
                        example, look in ml_tree/output for output trees for
                        ML tree inference)
        tree_names: list of tree files to look for in tree_base_path
        tree_fmt: ete3 Tree format to use for reading trees. Should be integers.
                  List, indiciating which tree format should be applied to read each
                  file in tree_names. Should be the same length as tree_names.
        df: topiary dataframe to use to interpret tip_columns. overrides the
            dataframe that is pulled from prev_run_dir.

    Return
    ------
        ete3 Tree instances for each tree requested in tree_names.
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
        df: topiary dataframe
        tip_columns: columns in dataframe to use to create human-readable names for
                     tips on tree. For example, tip_columns = ["species","nickname"]
                     would give names like "Homo sapiens|LY96". If None, try
                     ["species","nickname"]. If nickname is not in dataframe,
                     fall back to ["species",f"{name[:10]}..."]
        separator: separator to use between pretty names

    Return
    ------
        name_dict
    """

    local_df = df
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

    # Construct name_dict, which will map between the uid on the leaves to a
    # prettier/more useful name
    name_dict = {}
    for i in range(len(local_df)):
        key = local_df.uid.iloc[i]

        v = []
        for column in tip_columns:
            v.append(f"{local_df[column].iloc[i]}")
        value = separator.join(v)

        name_dict[key] = value

    return name_dict


def setup_generic_tree_format(num_leaves,
                              value_range=(0,100),
                              R=(1,0),G=(1,0),B=(1,1),
                              width=250,space_between_taxa=10,line_width=2):
    """
    Set up formating for a generic tree.

    Parameters
    ----------
        num_leaves: number of leaves in the tree
        value_range: value range for color map. if None, do not create a ColorMap.
        R: tuple holding red channel value corresponding to minimum and maximum
           values of span. Values should be between 0 and 1.
        G: tuple holding green channel value corresponding to minimum and maximum
           values of span. Values should be between 0 and 1.
        B: tuple holding blue channel value corresponding to minimum and maximum
           values of span. Values should be between 0 and 1.
        width: width of total tree (pixels)
        height: height of total tree (pixels)
        line_width: width of tree lines (pixels)

    Returns
    -------
        cm (color map)
        ts: ete3.TreeStyle object
        ns: ete3.NodeStyle object
    """

    # Decide whether/how we're doing color map
    if value_range is not None:
        cm = ColorMap(span=value_range,R=R,G=G,B=B)
    else:
        cm = None

    # Construct basic tree format
    ts = ete3.TreeStyle()
    ts.show_scale = True
    ts.scale = width
    ts.branch_vertical_margin = space_between_taxa

    # Construct basic node format
    ns = ete3.NodeStyle()
    ns["size"] = 0
    ns["hz_line_width"] = line_width
    ns["vt_line_width"] = line_width

    return cm, ts, ns

def final_render(T,ts,output_file,default_pdf):
    """
    Render tree T in a sane fashion.

    Parameters
    ----------
        T: ete3.Tree object to draw
        ts: ete3.TreeStyle object to use to format tree object
        output_file: write tree to output file. If None and in a notebook,
                     just draw on the screen. If None and *not* in a notebook,
                     write only to default_pdf
        default_pdf: write to this pdf file. if None, do not write.

    Returns
    -------
        T.render output (if in a notebook) or None (if not in a notebook. )
    """

    # If output_file is not specified...
    if output_file is None:

        # If not in a notebook, create a pdf
        if not topiary._in_notebook:
            output_file = default_pdf

    # If output_file specified, write out
    if output_file is not None:
        T.render(output_file,tree_style=ts)

    # If this is in a notebook, render it and return so it appears
    if topiary._in_notebook:
        return T.render("%%inline",tree_style=ts)

    return None
