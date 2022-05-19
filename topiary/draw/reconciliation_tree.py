__description__ = \
"""
Function for plotting tree from a species/gene tree reconciliation calculation.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-05-19"

import topiary
from topiary.external.interface import read_previous_run_dir

import os

from .base import load_trees, create_name_dict, setup_generic_tree_format, final_render
import ete3

def reconciliation_tree(run_dir,
                        output_file=None,
                        tip_columns=["species","nickname"],
                        tip_name_separator="|",
                        fontsize=36,
                        circle_radius=0.025,
                        event_color_map={"D":"green","L":"red","S":"blue","T":"pink"},
                        width=800,space_between_taxa=10,line_width=4,
                        df=None):
    """
    Draw a reconciliation tree, labeling the internal nodes according to
    whether they are speciations, duplications, or losses.

    Parameters
    ----------
        reconciliation_run_dir: directory containing a genrax reconcilation run
        df: dataframe (overides whatever is in ml_run_dir/output)
        output_file: output file to write tree. If running in a notebook, will
                     return to notebook and write to this file. If not in a notebook
                     and not specified, will write out reconciliation-tree.pdf.
        tip_columns: label the tree tips as "|".join(tip_columns). If tip_columns
                     is ["species","nickname"], tips will have names like
                     'Homo sapiens|LY96'.
        tip_name_separator: string to separate columns in tip names ("|" in
                            tip_columns example above.)
        fontsize: fontsize in points for tip labels and support values
        circle_radius: circle size for internal nodes (fraction of total width)
        event_color_map: map between events (duplication (D), loss (L), speciation (S),
                         or transfer (T) and node color).
        width: width of total tree (pixels)
        space_between_taxa: number of pixels between taxa, sets hight.
        line_width: width of tree lines (pixels)
        df: dataframe (overides whatever is in run_dir/output/dataframe.csv)

    Return
    ------
        None or rendered output if in notebook
    """

    # Load data from previous run
    prev_run = read_previous_run_dir(run_dir)

    # Load the tree and relevant information from the ml_run_dir output
    # directory.
    T = load_trees(prev_run_dir=run_dir,
                   tree_base_path=os.path.join("output","reconcilations"),
                   tree_names=["reconcile_events.newick"],
                   tree_fmt=[1])

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
                                           value_range=None,
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

    # Traverse T1 and T2 simultaneously
    for n in T.traverse():

        # Set base style on node of T1
        n.set_style(ns)

        # If this is not a leaf
        if not n.is_leaf():

            try:
                rgb = event_color_map[n.name]
            except KeyError:
                rgb = "gray"

            # Draw circles
            node_circle = ete3.CircleFace(radius=circle_radius,color=rgb,style="circle")
            node_circle.margin_right=-circle_radius
            n.add_face(node_circle,0,position="branch-right")

            if has_supports:

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

    return final_render(T,ts,output_file,"reconcilation-tree.pdf")
