"""
Draw a maximum likelihood tree with nodes colored by their bootstrap support.
"""

import topiary
from topiary.external._interface import read_previous_run_dir
from ._base import load_trees, create_name_dict, setup_generic_tree_format, final_render

import ete3

def ml_tree(run_dir,
            output_file=None,
            tip_columns=["species","nickname"],
            tip_name_separator="|",
            fontsize=36,
            circle_radius=0.025,
            value_range=(0,100),R=(0.95,0),G=(0.95,0),B=(1,1),
            width=800,space_between_taxa=10,line_width=4,
            df=None):
    """
    Draw a maximum likelihood tree with nodes colored by their bootstrap support.

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

    # Load the tree and relevant information from the run_dir output
    # directory.
    T = load_trees(prev_run_dir=run_dir,
                   tree_base_path="output",
                   tree_names=["tree.newick"],
                   tree_fmt=[0],
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
    cm, ts, ns = setup_generic_tree_format(num_leaves=len(T.get_leaves()),
                                           value_range=value_range,
                                           R=R,G=G,B=B,
                                           width=width,
                                           space_between_taxa=space_between_taxa,
                                           line_width=line_width)


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

    circle_radius = circle_radius*width

    # Traverse T
    for n in T.traverse():

        # Set base style on node of T
        n.set_style(ns)

        # If this is not a leaf
        if not n.is_leaf():

            if has_supports:

                if cm is not None:

                    # Get color for node
                    rgb = cm.hex(float(n.support))

                    # Draw circles on internal nodes
                    node_circle = ete3.CircleFace(radius=circle_radius,color=rgb,style="circle")
                    node_circle.margin_right=-circle_radius
                    n.add_face(node_circle,0,position="branch-right")

                # Wwrite out upport labels
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

    return final_render(T,ts,output_file,"ml-tree.pdf")
