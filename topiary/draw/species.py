"""
Draw a species tree with tips labeled by species name.
"""

import topiary
from .prettytree import PrettyTree

import ete3


def species_tree(species_tree,
                 output_file=None,
                 font_size=15,
                 stroke_width=2,
                 vertical_pixels_per_tip=20,
                 min_height=300,
                 tip_labels_align=True,
                 return_canvas=False,
                 **kwargs):
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
    font_size : float, default=15
        font size in points for labels
    stroke_width : int, default=2
        width of lines drawing tree (pixels)
    vertical_pixels_per_tip : int, default=20
        number of pixels to assign to each tip when calculating figure
        height
    min_height : float, default=300
        minimum height for figure (pixels)
    tip_labels_align : bool, default=True
        align species names on the right of the plot
    return_canvas : bool, default=False
        if True, return the canvas whether this is being called in a notebook
        or not. (canvas will be returned in a jupyter notebook regardless of this
        setting)
    **kwargs : dict, optional
        pass any other keyword arguments directly to toytree.tree.draw

    Returns
    -------
    plot : toyplot.canvas or None
        if running in jupyter notebook, return toyplot.canvas; otherwise, return
        None.
    """

    name_dict = {}
    for n in species_tree.traverse():
        if n.is_leaf():
            name_dict[n.name] = n.species

    pt = PrettyTree(species_tree,
                    name_dict=name_dict,
                    font_size=font_size,
                    stroke_width=stroke_width,
                    vertical_pixels_per_tip=vertical_pixels_per_tip,
                    min_height=min_height,
                    tip_labels_align=tip_labels_align,
                    **kwargs)

    # Figure out plotting to file.
    if output_file is not None:
        pt.render(output_file)

    # If this is in a notebook, return it so it appears
    ret = None
    if topiary._in_notebook or return_canvas:
        ret = pt.canvas

    return ret
