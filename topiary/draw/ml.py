"""
Draw a maximum likelihood tree with nodes colored by their bootstrap support.
"""

import topiary
from topiary._private.interface import read_previous_run_dir
from .core import load_trees, create_name_dict, final_render
from .prettytree import PrettyTree

import ete3

def ml_tree(run_dir,
            output_file=None,
            tip_columns=["species","nickname"],
            tip_name_separator="|",
            support_span=(50,100),
            color=("#ffffff","#155677"),
            size=None,
            font_size=15,
            stroke_width=2,
            vertical_pixels_per_taxon=20,
            min_height=300,
            df=None,
            **kwargs):
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
        'Homo sapiens|LY96'. If the name is not unique, the uid will be
        appended to the name.
    tip_name_separator : str, default="|"
        string to separate columns in tip names ("|" in tip_columns example
        above.) Cannot be "#,;:'\")(" as these are used in newick format.
    support_span : tuple, default=(50,100)
        set min/max values for support color calculation. First element
        is min, second is max. If either is None, take min or max from the
        min/max of the property
    color : str or tuple or dict, default=("#ffffff","#155677")
        set node color. If a single value, color all nodes that color. If
        list-like and length 2, treat as colors for minimum and maximum of a
        color gradient.  If dict, map property keys to color values. Colors
        can be RGBA tuples, named colors, or hexadecimal strings. See the
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

    # Make sure output file is a string
    if output_file is not None:
        output_file = str(output_file)

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

    # Rename leaves and grab supports
    supports = []
    for n in T.traverse():
        if not n.is_leaf():
            supports.append(n.support)

    # Figure out if the tree has supports. If the tree did not have supports,
    # ete3 will load it in with a support value of 1.0 for all nodes. This would
    # be very unlikely for a tree with real supports, so use that to look for
    # trees without supports
    supports_seen = list(set(supports))
    if len(supports_seen) == 1 and supports_seen[0] == 1:
        has_supports = False
    else:
        has_supports = True

    # Create tree
    pt = PrettyTree(T,
                    name_dict=name_dict,
                    font_size=font_size,
                    stroke_width=stroke_width,
                    vertical_pixels_per_taxon=vertical_pixels_per_taxon,
                    min_height=min_height,
                    **kwargs)
    pt.draw_scale_bar()

    # Draw supports
    if has_supports:
        pt.draw_nodes(property_label="support",
                      prop_span=support_span,
                      color=color,
                      size=size)
        pt.draw_node_labels(property_labels="support")
        pt.draw_node_legend()

    return final_render(pt,output_file=output_file,default_file="ml-tree.pdf")
