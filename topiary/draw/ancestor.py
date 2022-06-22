"""
Draw a tree with ancestors colored by their relative posterior probability.
"""

import topiary
from topiary.external._interface import read_previous_run_dir

import os, re, copy

from ._base import load_trees, create_name_dict, setup_generic_tree_format, final_render
import ete3

def ancestor_tree(run_dir,
                  output_file=None,
                  tip_columns=["species","nickname"],
                  tip_name_separator="|",
                  fontsize=36,
                  circle_radius=0.025,
                  value_range=(0.6,1),R=(0.95,0),G=(0.95,0),B=(1,1),
                  width=800,space_between_taxa=10,line_width=4,
                  df=None):
    """
    Draw a tree with ancestors colored by their relative posterior probability.

    Parameters
    ----------
    run_dir : str
        directory containing an ancestral reconstruction run
    output_file : str, optional
        output file to write tree. If running in a notebook, will return to
        notebook and write to this file. If running in a notebook but not
        specified, do not write to a file. If not in a notebook and not
        specified, will write out ancestor-tree.pdf.
    tip_columns: list, default=["species","nickname"]
        label the tree tips as "|".join(tip_columns). For example, if
        tip_columns is ["species","nickname"], tips will have names like
        'Homo sapiens|LY96'.
    tip_name_separator : str, default="|"
        string to separate columns in tip names ("|" in tip_columns example
        above.)
    fontsize : float, default=36
        fontsize in points for labels
    circle_radius : float, default=0.025
        circle size for internal nodes (fraction of total width)
    value_range : tuple,None default=(0.6,1)
        tuple holding minimum and maximum values for support color map.
        (0,100) means miminum support value is 0, maximum support value
        is 100. If None, disable node color map.
    R : tuple, default=(0.95,0)
        tuple holding red channel values for lowest and highest supports.
    G : tuple, default=(0.95,0)
        tuple holding green channel values for lowest and highest supports.
    B : tuple, default=(1.0,1.0)
        tuple holding blue channel values for lowest and highest supports.
    width : int, default=800
        width of total tree in pixels
    space_between_taxa : int, default=10
        number of pixels between taxa, sets height.
    line_width : int, default=4
        width of lines used to draw tree (pixels)
    df : pandas.DataFrame or None, default=None
        topiary dataframe (overides whatever is in run_dir/output/dataframe.csv)

    Returns
    -------
    Python.core.display.Image or None
        if running in jupyter notebook, return Image; otherwise, return None
    """

    # Load data from previous run
    prev_run = read_previous_run_dir(run_dir)

    # Load the tree and relevant information from the ml_run_dir output
    # directory.
    T_label, T_pp = load_trees(prev_run_dir=run_dir,
                               tree_base_path=os.path.join("output","ancestors"),
                               tree_names=["ancestors_label.newick",
                                           "ancestors_pp.newick"],
                               tree_fmt=[1,1],
                               outgroup=prev_run["outgroup"])

    # If df not specified, get from the previous run
    if df is None:
        df = prev_run["df"]

    # Create dictionary mapping between uid and pretty name format
    name_dict = create_name_dict(df=df,
                                 tip_columns=tip_columns,
                                 separator=tip_name_separator)

    # Set up formats (color map for main nodes, tree format (ts) and generic
    # node format (ns)).
    cm, ts, ns = setup_generic_tree_format(num_leaves=len(T_label.get_leaves()),
                                           value_range=value_range,
                                           R=R,G=G,B=B,
                                           width=width,
                                           space_between_taxa=space_between_taxa,
                                           line_width=line_width)


    circle_radius = circle_radius*width

    # Grab the root from these trees. If the root as a posterior probability,
    # the annotation shifted when we rooted the tree. Grab the annotation and
    # put it on the inappropriately empty node that should have this annoation
    label_root = T_label.get_tree_root()
    pp_root = T_pp.get_tree_root()
    if pp_root.name != "":
        offset_root = (copy.deepcopy(label_root.name),copy.deepcopy(pp_root.name))
    else:
        offset_root = None

    # Traverse posterior probability and label trees simultaneously
    T_label_traversed = list(T_label.traverse())
    T_pp_traversed = list(T_pp.traverse())
    for i in range(len(T_label_traversed)):

        # Get T1 and T2 nodes (n and m)
        n = T_label_traversed[i]
        m = T_pp_traversed[i]

        # Set base style on node of T1
        n.set_style(ns)

        # Do not draw anything for root on ancestral reconstruction tree
        if n.is_root():
            continue

        # If this is not a leaf
        if not n.is_leaf():

            anc_name = n.name
            pp = m.name

            # If posterior probability is not defined for a node that is *not*
            # the root...
            if pp == "":

                # If the root was mislabeled above, stick the label on this node
                if offset_root is not None:
                    anc_name = offset_root[0]
                    pp = offset_root[1]
                    offset_root = None

            if cm is not None:

                try:
                    v = float(pp)
                    rgb = cm.hex(v)

                    # Draw circles
                    node_circle = ete3.CircleFace(radius=circle_radius,color=rgb,style="circle")
                    node_circle.margin_right=-circle_radius
                    n.add_face(node_circle,0,position="branch-right")

                except ValueError:
                    pass

            # Ancestor labels
            anc_label = ete3.TextFace(text=re.sub("anc","a",anc_name),fsize=fontsize)
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
            txt = ete3.TextFace(clean_name,fsize=fontsize)
            txt.margin_left = 4
            n.add_face(txt,0,position="branch-right")


    return final_render(T_label,ts,output_file,"ancestor-tree.pdf")
