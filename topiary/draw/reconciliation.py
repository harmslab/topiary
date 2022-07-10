"""
Draw a reconciliation tree, labeling the internal nodes according to
whether they are speciations, duplications, or losses.
"""

import topiary
from topiary.external._interface import read_previous_run_dir

import os

from .core import load_trees, create_name_dict, final_render
from .prettytree import PrettyTree
import ete3

def reconciliation_tree(run_dir,
                        output_file=None,
                        tip_columns=["species","nickname"],
                        tip_name_separator="|",
                        event_color_map={"D":"green","L":"red","S":"blue","T":"pink"},
                        font_size=15,
                        stroke_width=2,
                        vertical_pixels_per_taxon=20,
                        min_height=300,
                        df=None,
                        **kwargs):
    """
    Draw a reconciliation tree, labeling the internal nodes according to
    whether they are speciations, duplications, or losses.

    Parameters
    ----------
    run_dir : str
        directory containing an ancestral reconstruction run
    output_file : str, optional
        output file to write tree. If running in a notebook, will return to
        notebook and write to this file. If running in a notebook but not
        specified, do not write to a file. If not in a notebook and not
        specified, will write out ancestor-tree.pdf.
    tip_columns : list, default=["species","nickname"]
        label the tree tips as "|".join(tip_columns). For example, if
        tip_columns is ["species"|"nickname"], tips will have names like
        'Homo sapiens|LY96'. If the name is not unique, the uid will be
        appended to the name.
    tip_name_separator : str, default="|"
        string to separate columns in tip names ("|" in tip_columns example
        above.) Cannot be "#,;:'\")(" as these are used in newick format.
    event_color_map : dict
        map between evolutionary events and colors. Keys should be D,L,S, and T.
        color values can be RGBA, named colors, or hexadecimal strings. See the
        toyplot documentation for allowed values.
    font_size : float, default=15
        font size in points for labels
    stroke_width : int, default=2
        width of lines drawing tree (pixels)
    vertical_pixels_per_taxon : int, default=20
        number of pixels to assign to each taxon when calculating figure
        height
    min_height : float, default=300
        minimum height for figure (pixels)
    df : pandas.DataFrame or None, default=None
        topiary dataframe (overides whatever is in run_dir/output/dataframe.csv)
    **kwargs : dict, optional
        pass any other keyword arguments directly to toytree.tree.draw

    Returns
    -------
    plot : toyplot.canvas or None
        if running in jupyter notebook, return toyplot.canvas; otherwise, return
        None.
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

    # Give pretty names to tips
    for n in T.traverse():
        if n.is_leaf():
            n.name = name_dict[n.name]

    # Create tree
    pt = PrettyTree(T,
                    font_size=font_size,
                    stroke_width=stroke_width,
                    vertical_pixels_per_taxon=vertical_pixels_per_taxon,
                    min_height=min_height,
                    **kwargs)
    pt.draw_scale_bar()
    pt.draw_nodes(property_label="name",
                  color=event_color_map)
    pt.draw_node_legend()

    return final_render(pt,output_file=output_file,default_file="reconcilation-tree.pdf")
