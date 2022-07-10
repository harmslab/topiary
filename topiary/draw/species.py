"""
Draw a species tree with tips labeled by species name.
"""

import topiary
from .core import final_render

import ete3


def species_tree(species_tree,
                 output_file=None,
                 fontsize=36,
                 circle_radius=0.0,
                 width=800,space_between_taxa=10,line_width=4):
    """
    Draw a species tree with tips labeled by species name.

    Parameters
    ----------
    species_tree: ete3.Tree
        ete3 tree object holding a species tree
    output_file : str, optional
        output file to write tree. If running in a notebook, will return to
        notebook and write to this file. If running in a notebook but not
        specified, do not write to a file. If not in a notebook and not
        specified, will write out ancestor-tree.pdf.
    fontsize : float, default=36
        fontsize in points for labels
    circle_radius : float, default=0.025
        circle size for internal nodes (fraction of total width)
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

    T = species_tree.copy()

    # Set up formats (color map for main nodes, tree format (ts) and generic
    # node format (ns)).
    cm, ts, ns = setup_generic_tree_format(num_leaves=len(T.get_leaves()),
                                           width=width,
                                           space_between_taxa=space_between_taxa,
                                           line_width=line_width)

    for n in T.traverse():

        n.set_style(ns)

        if n.is_leaf():

            # Add text for clean name
            n.name = ""
            txt = ete3.TextFace(n.species[0],fsize=fontsize)
            txt.margin_left = 4
            n.add_face(txt,0,position="branch-right")


    return final_render(T,ts,output_file,"species-tree.pdf")
