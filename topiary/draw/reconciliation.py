"""
Draw a reconciliation tree, labeling the internal nodes according to
whether they are speciations, duplications, or losses.
"""

import topiary

from topiary._private.interface import read_previous_run_dir
from ._core import load_trees, create_name_dict, final_render, map_tree_to_tree
from ._core import draw_bs_nodes
from ._prettytree import PrettyTree

import numpy as np

import os


def reconciliation_tree(run_dir,
                        output_file=None,
                        event_color={"D":"#64007F","L":"#BAD316","S":"#D16A16","T":"#407E98"},
                        bs_color={50:"#ffffff",100:"#155677"},
                        event_label=False,
                        bs_label=False,
                        tip_columns=["species","nickname"],
                        tip_name_separator="|",
                        color=None,
                        size=None,
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
    size : float or tuple or dict, optional
        set node size. If a single value, make all nodes that size. If
        list-like and length 2, treat as sizes for minimum and maximum of a
        size gradient. If dict, map property keys to size values. Sizes must
        be float >= 0.
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

    # Figure out what kind of parsing/drawing to do
    if prev_run["calc_type"] == "reconciliation":
        plot_supports = False
    elif prev_run["calc_type"] == "reconciliation_bootstrap":
        plot_supports = True
    else:
        err = f"\nrun_dir '{run_dir}' is not a reconcilation calculation.\n\n"
        raise ValueError(err)

    # Load the events tree
    T_event = load_trees(prev_run_dir=run_dir,
                         tree_base_path=os.path.join("output","reconcilations"),
                         tree_names=["reconcile_events.newick"],
                         tree_fmt=[1])

    # Load tree with supports if there are supports in calculation
    T_support = None
    if plot_supports:
        T_support = load_trees(prev_run_dir=run_dir,
                               tree_base_path="output",
                               tree_names=["tree.newick"],
                               tree_fmt=[1])

        shared, T_event_only, T_support_only= map_tree_to_tree(T_event,
                                                               T_support)
    else:
        shared = []
        T_event_only = list(T_event.traverse())
        T_support_only = []

    # If df not specified, get from the previous run
    if df is None:
        df = prev_run["df"]

    # Create dictionary mapping between uid and pretty name format
    name_dict = create_name_dict(df=df,
                                 tip_columns=tip_columns,
                                 separator=tip_name_separator)

    # Go through shared nodes and label T_event nodes with support. (Shared will
    # have len(shared) == 0 if no supports loaded)
    for i in range(len(shared)):
        n = shared[i][0]
        if not n.is_leaf():
            m = shared[i][1]
            if m.name != "":
                n.support = float(m.name)
            n.add_feature("event",n.name)

    # Create tree
    pt = PrettyTree(T_event,
                    name_dict=name_dict,
                    font_size=font_size,
                    stroke_width=stroke_width,
                    vertical_pixels_per_taxon=vertical_pixels_per_taxon,
                    min_height=min_height,
                    **kwargs)
    pt.draw_scale_bar()

    added_supports = False
    if plot_supports:
        added_supports, size = draw_bs_nodes(T_event,
                                             pt,
                                             bs_color,
                                             bs_label,
                                             color,
                                             size)

    # color takes precedence over event_color.
    if color is not None:
        event_color = color

    # if event_color is defined (either by "color" or by event_color argument)
    if event_color is not None:
        pt.draw_nodes(property_label="event",
                      color=event_color,
                      size=size)

    pt.draw_node_legend()

    return final_render(pt,output_file=output_file,default_file="reconcilation-tree.pdf")
