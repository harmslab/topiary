__description__ = \
"""
Functions for plotting ete3 trees.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-05-12"

import topiary
import ete3

import numpy as np
import re, os

class ColorMap:

    def __init__(self,span=(0,1),R=(1,0),G=(1,0),B=(1,0)):

        self._min_value = span[0]
        self._max_value = span[1]

        self._R_bottom = R[0]
        self._R_top = R[1]

        self._R_slope = (self._R_top - self._R_bottom)
        self._R_intercept = self._R_bottom

        self._G_bottom = G[0]
        self._G_top = G[1]
        self._G_slope = (self._G_top - self._G_bottom)
        self._G_intercept = self._G_bottom

        self._B_bottom = B[0]
        self._B_top = B[1]
        self._B_slope = (self._B_top - self._B_bottom)
        self._B_intercept = self._B_bottom

    def hex(self,value):

        if value < self._min_value:
            value = self._min_value
        if value > self._max_value:
            value = self._max_value

        scaled = (value - self._min_value)/(self._max_value - self._min_value)

        R = int(round(self._R_intercept + self._R_slope*scaled,0))
        G = int(round(self._G_intercept + self._G_slope*scaled,0))
        B = int(round(self._B_intercept + self._B_slope*scaled,0))

        return f"#{bytearray((R,G,B)).hex()}"

def _load_trees(prev_dir,
                tree_base_path,
                tree_names,
                tip_columns,
                tree_fmt=None,
                df=None):
    """
    """

    # Load previous run
    prev_info = topiary.io.read_previous_run_dir(prev_dir)

    # If user did not pass in a dataframe, grab dataframe from previous run
    if df is None:
        df = prev_info["df"]

    # Make sure the tip columns are sane
    for c in tip_columns:
        try:
            df[c]
        except KeyError:
            err = f"\ncolumn '{c}' not in dataframe. Please check the\n"
            err += "tip_columns argument.\n\n"
            raise ValueError(err)

    # Construct name_dict, which will map between the uid on the leaves to a
    # prettier/more useful name
    name_dict = {}
    for i in range(len(df)):
        key = df.uid.iloc[i]

        v = []
        for column in tip_columns:
            v.append(f"{df[column].iloc[i]}")
        value = "|".join(v)

        name_dict[key] = value

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

        tree_path = os.path.join(tree_base_path,t)
        if not os.path.exists(tree_path):
            err = f"\nCould not find tree file '{tree_path}'\n\n"
            raise FileNotFoundError(err)

        T_list.append(topiary.io.read_tree(tree_path,fmt=tree_fmt[i]))

    return name_dict, *T_list


def _setup_generic_tree_format(value_range=(0,100),
                               R=(255,0),G=(255,0),B=(255,0),
                               x_size=250,y_size=200,line_width=2):

    # Decide whether/how we're doing color map
    if value_range is not None:
        cm = ColorMap(span=value_range,R=R,G=G,B=B)
    else:
        cm = None

    # Construct basic tree format
    ts = ete3.TreeStyle()
    ts.show_scale = True
    ts.scale = x_size
    ts.branch_vertical_margin = y_size/len(T.get_leaves())

    # Construct basic node format
    ns = ete3.NodeStyle()
    ns["size"] = 0
    ns["hz_line_width"] = line_width
    ns["vt_line_width"] = line_width

    return cm, ts, ns

def _final_render(T,ts,output,default_pdf):

    # If output is not specified...
    if output is None:

        # If not in a notebook, create a pdf
        if not topiary._in_notebook:
            output = default_pdf

    # If output specified, write out
    if output is not None:
        T.render(output,tree_style=ts)

    # If this is in a notebook, render it and return so it appears
    if topiary._in_notebook:
        return T.render("%%inline",tree_style=ts)

    return None

def species_tree(species_tree,output=None):
    """
    species_tree: ete3 Tree object holding a species tree
    output: None or output file. If None and work being done
            in a notebook, render inline. If not in a notebook,
            write to species_tree.pdf.
    """

    T = species_tree.copy()

    ts = ete3.TreeStyle()
    ts.draw_guiding_lines = True
    ts.guiding_lines_type = 0
    ts.show_scale = False
    ts.scale = 20

    for n in T.traverse():

        ns = ete3.NodeStyle()
        ns["size"] = 0
        ns["hz_line_width"] = 2
        ns["vt_line_width"] = 2

        n.set_style(ns)

        if n.is_leaf():
            n.name = ""

            txt = ete3.TextFace(n.species[0],fsize=14)
            txt.margin_left = 4
            n.add_face(txt,0,position="branch-right")


    if output is None:
        if not topiary._in_notebook:
            output = "species_tree.pdf"

    if output is not None:
        T.render(output,tree_style=ts)

    if topiary._in_notebook:
        return T.render("%%inline",tree_style=ts)

def ml_tree(ml_dir,
            df=None,
            output=None,
            tip_columns=["species","nickname"],
            fontsize=14,
            circle_radius=5,
            value_range=(0,100),R=(255,0),G=(255,0),B=(255,0),
            x_size=250,y_size=200,line_width=2):
    """
    """

    name_dict, T = _load_trees(prev_dir="ml_dir",
                               tree_base_path="output",
                               tree_names=["tree.newick"],
                               tree_fmt=[0],
                               tip_columns=tip_columns,
                               df=df)

    cm, ts, ns = _setup_generic_tree_format(value_range=value_range,
                                            R=R,G=G,B=B,
                                            x_size=x_size,y_size=y_size,
                                            line_width=lnie_width)


    # Figure out if the tree has supports. If the tree did not have supports,
    # ete3 will load it in with a support value of 1.0 for all nodes. This would
    # be very unlikely for a tree with real supports, so use that to look for
    # trees without supports
    supports = []
    for n in T.traverse():
        if not n.is_leaf():
            supports.append(n.support)

    supports_seen = list(set(supports))
    if len(supports_seen) == 1 and supports_seen[0] == 1:
        has_supports = False
    else:
        has_supports = True

    # Traverse Ts
    for n in T.traverse():

        # Set base style on node of T
        n.set_style(ns)

        # If this is not a leaf
        if not n.is_leaf():

            if has_supports:

                if cm is not None:

                    rgb = cm.hex(float(n.support))

                    # Draw circles
                    node_circle = ete3.CircleFace(radius=circle_radius,color=rgb,style="circle")
                    node_circle.margin_right=-circle_radius
                    n.add_face(node_circle,0,position="branch-right")

                # Support labels
                anc_label = ete3.TextFace(text=f"{int(round(n.support,0)):d}",fsize=fontsize)
                anc_label.inner_background.color = "white"
                anc_label.margin_right = circle_radius
                anc_label.margin_top = 0
                n.add_face(anc_label,0,position="float")

        # If this is a leaf
        else:

            # Get the clean name (rather than uid)
            clean_name = name_dict[n.name]

            # Add text for clean name
            n.name = ""
            txt = ete3.TextFace(clean_name,fsize=fontsize)
            txt.margin_left = 4
            n.add_face(txt,0,position="branch-right")


    return _final_render((T,ts,output,"ml-tree.pdf")

def reconciliation_tree(reconcilation_dir,
                        df=None,
                        output=None,
                        tip_columns=["species","nickname"]):
    """
    """

    # Load previous run
    rec_info = topiary.io.read_previous_run_dir(reconcilation_dir)

    # If user did not pass in a dataframe, grab dataframe from previous run
    if df is None:
        df = rec_info["df"]

    # Load the output trees from the calculation
    tree_path = os.path.join(reconcilation_dir,"output","reconcilations")
    T_path = os.path.join(tree_path,"reconcile_events.newick")

    if not os.path.exists(T_path):
        err = f"\nreconcilation directory '{reconcilation_dir}' does not have\n"
        err += "reconcile_events.newick.\n"
        raise FileNotFoundError(err)

    T = topiary.io.read_tree(T_path,fmt=1)

    # Make sure the tip columns are sane
    for c in tip_columns:
        try:
            df[c]
        except KeyError:
            err = f"\ncolumn '{c}' not in dataframe. Please check the\n"
            err += "tip_columns argument.\n\n"
            raise ValueError(err)

    # Construct name_dict, which will map between the uid on the leaves to a
    # prettier/more useful name
    name_dict = {}
    for i in range(len(df)):
        key = df.uid.iloc[i]

        v = []
        for column in tip_columns:
            v.append(f"{df[column].iloc[i]}")
        value = "|".join(v)

        name_dict[key] = value

    # Construct basic tree format
    ts = ete3.TreeStyle()
    ts.show_scale = True
    ts.scale = 250
    ts.branch_vertical_margin = 5

    # Construct basic node format
    ns = ete3.NodeStyle()
    ns["size"] = 0
    ns["hz_line_width"] = 2
    ns["vt_line_width"] = 2

    # Traverse T1 and T2 simultaneously
    for n in T.traverse():

        # Set base style on node of T1
        n.set_style(ns)

        # If this is not a leaf
        if not n.is_leaf():

            if n.name == "D":
                rgb = "red"
            else:
                rgb = "blue"

            # Draw circles
            node_circle = ete3.CircleFace(radius=5,color=rgb,style="circle")
            node_circle.margin_right=-5
            n.add_face(node_circle,0,position="branch-right")

        # If this is a leaf
        else:

            # Get the clean name (rather than uid)
            clean_name = name_dict[n.name]

            # Add text for clean name
            n.name = ""
            txt = ete3.TextFace(clean_name,fsize=14)
            txt.margin_left = 4
            n.add_face(txt,0,position="branch-right")


    # If output is not specified...
    if output is None:

        # If not in a notebook, create a pdf
        if not topiary._in_notebook:
            output = "reconcilation_tree.pdf"

    # If output specified, write out
    if output is not None:
        T.render(output,tree_style=ts)

    # If this is in a notebook, render it and return so it appears
    if topiary._in_notebook:
        return T.render("%%inline",tree_style=ts)

def ancestor_tree(ancestor_dir,
                  df=None,
                  output=None,
                  tip_columns=["species","nickname"]):
    """
    """

    # Load previous run
    anc_info = topiary.io.read_previous_run_dir(ancestor_dir)

    # If user did not pass in a dataframe, grab dataframe from previous run
    if df is None:
        df = anc_info["df"]

    # Load the output trees from the calculation
    tree_path = os.path.join(ancestor_dir,"output","ancestors")
    T1_path = os.path.join(tree_path,"ancestors_label.newick")
    T2_path = os.path.join(tree_path,"ancestors_pp.newick")

    if not os.path.exists(T1_path) or not os.path.exists(T2_path):
        err = f"\nancestors directory '{ancestor_dir}' does not have\n"
        err += "ancestors_label.newick or ancestors_pp.newick.\n"
        raise FileNotFoundError(err)

    T1 = topiary.io.read_tree(T1_path,fmt=1)
    T2 = topiary.io.read_tree(T2_path,fmt=1)

    # Make sure the tip columns are sane
    for c in tip_columns:
        try:
            df[c]
        except KeyError:
            err = f"\ncolumn '{c}' not in dataframe. Please check the\n"
            err += "tip_columns argument.\n\n"
            raise ValueError(err)

    # Construct name_dict, which will map between the uid on the leaves to a
    # prettier/more useful name
    name_dict = {}
    for i in range(len(df)):
        key = df.uid.iloc[i]

        v = []
        for column in tip_columns:
            v.append(f"{df[column].iloc[i]}")
        value = "|".join(v)

        name_dict[key] = value

    # Construct basic tree format
    ts = ete3.TreeStyle()
    ts.show_scale = True
    ts.scale = 250
    ts.branch_vertical_margin = 5

    # Construct basic node format
    ns = ete3.NodeStyle()
    ns["size"] = 0
    ns["hz_line_width"] = 2
    ns["vt_line_width"] = 2

    # Traverse T1 and T2 simultaneously
    T1_traversed = list(T1.traverse())
    T2_traversed = list(T2.traverse())
    for i in range(len(T1_traversed)):

        # Get T1 and T2 nodes (n and m)
        n = T1_traversed[i]
        m = T2_traversed[i]

        # Set base style on node of T1
        n.set_style(ns)

        # If this is not a leaf
        if not n.is_leaf():

            # Create ancestor circles from white (0.6) to blue (1.0)
            pp = (float(m.name) - 0.6)/.4
            channel = int(np.round(255*pp,0))
            rgb = f"#{bytearray((255-channel,255-channel,channel)).hex()}"

            # Draw circles
            node_circle = ete3.CircleFace(radius=5,color=rgb,style="circle")
            node_circle.margin_right=-5
            n.add_face(node_circle,0,position="branch-right")

            # Ancestor labels
            anc_label = ete3.TextFace(text=re.sub("anc","a",n.name),fsize=14)
            anc_label.inner_background.color = "white"
            anc_label.margin_right = 5
            anc_label.margin_top = 0
            n.add_face(anc_label,0,position="float")

        # If this is a leaf
        else:

            # Get the clean name (rather than uid)
            clean_name = name_dict[n.name]

            # Add text for clean name
            n.name = ""
            txt = ete3.TextFace(clean_name,fsize=14)
            txt.margin_left = 4
            n.add_face(txt,0,position="branch-right")


    # If output is not specified...
    if output is None:

        # If not in a notebook, create a pdf
        if not topiary._in_notebook:
            output = "ancestor_tree.pdf"

    # If output specified, write out
    if output is not None:
        T1.render(output,tree_style=ts)

    # If this is in a notebook, render it and return so it appears
    if topiary._in_notebook:
        return T1.render("%%inline",tree_style=ts)
