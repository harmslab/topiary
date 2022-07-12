"""
Draw a tree with ancestors colored by their relative posterior probability.
"""

import topiary
from topiary._private.interface import read_previous_run_dir

import os, re, copy

from ._core import load_trees, create_name_dict, final_render
from ._prettytree import PrettyTree
import ete3


def ancestor_tree(run_dir,
                  output_file=None,
                  tip_columns=["species","nickname"],
                  tip_name_separator="|",
                  support_span=(0.5,1.0),
                  color=("#ffffff","#D16A16"),
                  size=None,
                  font_size=15,
                  stroke_width=2,
                  vertical_pixels_per_taxon=20,
                  min_height=300,
                  df=None,
                  **kwargs):
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
        'Homo sapiens|LY96'. If the name is not unique, the uid will be
        appended to the name.
    tip_name_separator : str, default="|"
        string to separate columns in tip names ("|" in tip_columns example
        above.) Cannot be "#,;:'\")(" as these are used in newick format.
    support_span : tuple, default=(0.5,1.0)
        set min/max values for posterior probability calculation. First element
        is min, second is max. If None, take full span.
    color : str or tuple or dict, default=("#ffffff","#D16A16")
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

    # Grab the root from these trees. If the root has a posterior probability,
    # the annotation shifted when we rooted the tree. Grab the annotation and
    # put it on the inappropriately empty node that should have this annoation
    label_root = T_label.get_tree_root()
    pp_root = T_pp.get_tree_root()
    if pp_root.name != "":
        offset_root = (copy.deepcopy(label_root.name),copy.deepcopy(pp_root.name))
        pp_root.name = -1
        label_root.name = ""
    else:
        offset_root = None

    # Traverse posterior probability and label trees simultaneously
    T_label_nodes = list(T_label.traverse())
    T_pp_nodes = list(T_pp.traverse())
    for i in range(len(T_label_nodes)):

        # Get T1 and T2 nodes (n and m)
        n = T_label_nodes[i]
        m = T_pp_nodes[i]

        # If this is not a leaf
        if not n.is_leaf():

            pp = m.name

            # pp == "" only happens for root node
            if pp == "":

                # If the root was mislabeled above, stick the label on this node
                if offset_root is not None:
                    n.add_feature("name",offset_root[0])
                    n.add_feature("pp",offset_root[1])
                    offset_root = None

            else:
                n.add_feature("pp",m.name)


            n.name = re.sub("anc","a",n.name)

    pt = PrettyTree(T_label,
                    name_dict=name_dict,
                    font_size=font_size,
                    stroke_width=stroke_width,
                    vertical_pixels_per_taxon=vertical_pixels_per_taxon,
                    min_height=min_height,
                    **kwargs)
    pt.draw_scale_bar()

    # Draw supports
    pt.draw_nodes(property_label="pp",
                  prop_span=support_span,
                  plot_root=False,
                  color=color,
                  size=size)
    pt.draw_node_labels(property_labels="name")
    pt.draw_node_legend(label_renamer={"pp":"post. prob."})

    return final_render(pt,output_file=output_file,default_file="ancestor-tree.pdf")
