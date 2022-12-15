"""
Draw a topiary tree with nodes colored by calculation outputs.
"""

import topiary

from topiary._private import Supervisor
from topiary.draw.core import create_name_dict
from topiary.draw.core import construct_sizemap
from topiary.io.tree import load_trees
from topiary.draw.prettytree import PrettyTree

import toyplot
import numpy as np

import os, copy

def tree(calculation,
         output_file=None,
         bs_color={50:"#ffffff",100:"#000000"},
         pp_color={0.7:"#ffffff",1.0:"#DC801A"},
         event_color={"D":"#64007F","L":"#BAD316","T":"#407E98"}, #"S":"#023E55"}, #D16A16"},
         bs_label=False,
         pp_label=False,
         event_label=False,
         anc_label=True,
         tip_columns=None,
         tip_name_separator="|",
         disambiguate_tip_names=True,
         node_color=None,
         node_size=None,
         font_size=15,
         tip_text_style=None,
         label_text_style=None,
         label_position="right",
         label_position_offset=None,
         label_color="gray",
         stroke_width=2,
         vertical_pixels_per_tip=25,
         min_height=300,
         df=None,
         anc_link_path=None,
         return_canvas=False,
         **kwargs):
    """
    Draw a tree with annotated with calculation outputs.

    Parameters
    ----------
    calculation : str or Supervisor
        completed calculation. Should either be a directory containing the calc
        (e.g. the directory with run_parameters.json, input, working, output) or
        a Supervisor instance with a calculation loaded.
    output_file : str, optional
        output file to write tree.
    bs_color : dict, default={50:"#ffffff",100:"#000000"}
        set min/max values for branch support color map. First key is min,
        second is max. Values are valid toyplot colors (see Notes).
        If given, color argument takes precedence over bs_color.
    pp_color : dict, default={0.7:"#ffffff",1.0:"#DC801A"}
        set min/max values for posterior probability color map. First key is
        min, second is max. Values are valid toyplot colors (see Notes).
        If given, node_color argument takes precedence over pp_color.
    event_color : dict, default={"D":"#64007F","L":"#BAD316","T":"#407E98"}
        colors for evolutionary events determined by gene/species tree
        reconciliation. Allowed keys are "S" (speciation), "D" (duplication),
        "L" (loss), and "T" (transfer). If a key is not specified, that
        event will not be drawn on the tree. (The default is to *not* show
        speciation). Values are any valid toyplot colors (see Notes).
        If given, node_color argument takes precedence over event_color.
    bs_label : bool, default=False
        whether or not to write branch support values on the tree
    pp_label : bool, default=False
        whether or not to write posterior probability values on the tree
    event_label : bool, default=False
        whether or not to write events as text on the tree
    anc_label : bool, default=True
        whether or not to write ancestor names as text on the tree
    tip_columns: list, optional
        label the tree tips as "|".join(tip_columns). For example, if
        tip_columns is ["species","recip_paralog"], tips will have names like
        'Homo sapiens|LY96'. If the name is not unique, the uid will be
        appended to the name (see disambiguate_tip_names below). If not
        specified, try ["species","recip_paralog"], ["species","nickname"],
        then, finally, ["species","name"].
    tip_name_separator : str, default="|"
        string to separate columns in tip names ("|" in tip_columns example
        above.)
    disambiguate_tip_names : bool, default=True
        if two tip labels will be the same (for example, two labels will be
        "Homo sapiens|LY96"), append the uid to those labels so they can be
        uniquely identified.
    node_color : str or tuple or dict, optional
        internal node color. If a single value, color all nodes that color. If
        list-like and length 2, treat as colors for minimum and maximum of a
        color gradient.  If dict, map property keys to color values. This
        will override bs_color, pp_color, and event_color. We recommend
        using this to specify single colors only. Use bs_color, pp_color,
        and event_color to color nodes based on their properties.
    node_size : float or tuple or dict, optional
        set node size in pixels. If a single value, make all nodes that size.
        If list-like and length 2, treat as sizes for minimum and maximum of a
        size gradient. If dict, map property keys to size values. Sizes must
        be float >= 0. The size map will be applied for each property
        individually, so this could lead to strange results if plotting
        branch supports and posterior probabilities simultaneously.
    font_size : float, default=15
        font size in points for labels
    tip_text_style : dict, optional
        dictionary of css key/value pairs for drawing tip labels. (See Notes
        for details.)
    label_text_style : dict, optional
        dictionary of css key/value pairs for drawing internal nodes. (See
        Notes for details).
    label_position : str, default="right"
        where to put internal labels relative to node. allowed values are
        "top-left", "top", "top-right", "right", "bottom-right", "bottom",
        "bottom-left", and "left".
    label_position_offset : float, optional
        how far to displace the internal labels off the nodes in pixels.
    stroke_width : int, default=2
        width of lines drawing tree in pixels
    vertical_pixels_per_tip : int, default=20
        number of pixels to assign to each tip when calculating figure
        height
    min_height : float, default=300
        minimum height for figure (pixels)
    df : pandas.DataFrame or None, default=None
        topiary dataframe (overrides whatever is in run_directory/output/dataframe.csv)
    anc_link_path : str, optional
        if specified, format ancestors as links to ancestors in anc_link_path. 
    return_canvas : bool, default=False
        if True, return the canvas whether this is being called in a notebook
        or not. (canvas will be returned in a jupyter notebook regardless of this
        setting)
    **kwargs : dict, optional
        pass any other keyword arguments directly to toytree.tree.draw

    Notes
    -----
    + Colors can be specified as RGBA tuples, named colors, or hexadecimal strings.
      See the `toyplot documentation <https://toyplot.readthedocs.io/en/stable/colors.html>`_
      for allowed values.
    + Allowed css keys are determined by toyplot. See the
      `here <https://toyplot.readthedocs.io/en/stable/_modules/toyplot/style.html>`_
      for details.

    Returns
    -------
    plot : toyplot.canvas or None
        if running in jupyter notebook, return toyplot.canvas; otherwise, return
        None.
    """

    if issubclass(type(calculation),str):
        supervisor = Supervisor(calculation)
    else:
        supervisor = calculation

    # Make sure output file is a string
    if output_file is not None:
        output_file = str(output_file)

    # Load a tree with current calculation states (will have event, bs_support,
    # anc_label, anc_pp) on internal nodes. None if those parameters were not
    # calculated in at this point in the pipeline.
    T = load_trees(directory=supervisor.output_dir,
                   prefix=supervisor.tree_prefix)

    # If df not specified, get from the previous run
    if df is None:
        df = supervisor.df

    if tip_columns is None:
        if "recip_paralog" in df.columns:
            recip_columns = ["species","recip_paralog"]
        elif "nickname" in df.columns:
            recip_columns = ["species","nickname"]
        else:
            recip_columns = ["species","name"]

    # Create dictionary mapping between uid and pretty name format
    name_dict = create_name_dict(df=df,
                                 tip_columns=tip_columns,
                                 separator=tip_name_separator,
                                 disambiguate=disambiguate_tip_names)

    # Create tree
    pt = PrettyTree(T,
                    name_dict=name_dict,
                    font_size=font_size,
                    stroke_width=stroke_width,
                    vertical_pixels_per_tip=vertical_pixels_per_tip,
                    min_height=min_height,
                    **kwargs)

    # Set node_size if not set
    if node_size is None:
        node_size = pt.default_size

    # node_color takes precendence over specific entries
    if node_color is not None:
        pt.draw_nodes(color=node_color,size=node_size)

    else:

        # Plot events
        if event_color is not None:

            # Make sure node_size, whatever it is, works fine with event
            prop = []
            for n in T.traverse():
                if not n.is_leaf():
                    prop.append(n.event)
            prop = list(set(prop))
            sm, sm_span = construct_sizemap(node_size,prop)

            # Now update the size so it's slightly bigger than requested for
            # the event
            this_size = copy.deepcopy(node_size)
            if issubclass(type(this_size),dict):
                for k in this_size:
                    this_size[k] = this_size[k]*1.5
            else:
                try:
                    this_size = this_size*1.5
                except (ValueError,TypeError):
                    this_size = np.array(this_size)*1.5

            pt.draw_nodes(property_label="event",
                          color=event_color,
                          size=this_size)


        # Plot bootstrap supports
        if bs_color is not None:

            plot_bs, bs_span, bs_color = topiary.draw.core.parse_span_color(bs_color,node_color)
            pt.draw_nodes(property_label="bs_support",
                          prop_span=bs_span,
                          color=bs_color,
                          size=node_size)

            # If we successfully plotted bootstraps decrease node_size by factor
            # so we can plot next data
            if "bs_support" in pt.plotted_properties:

                if issubclass(type(node_size),dict):
                    for k in node_size:
                        node_size[k] = node_size[k]*0.6
                else:
                    try:
                        node_size = node_size*0.6
                    except (ValueError,TypeError):
                        node_size = np.array(node_size)*0.6

        # Plot ancestor posterior probabilities
        if pp_color is not None:
            plot_pp, pp_span, pp_color = topiary.draw.core.parse_span_color(pp_color,node_color)
            pt.draw_nodes(property_label="anc_pp",
                          prop_span=pp_span,
                          color=pp_color,
                          size=node_size)

    # Deal with property labeling
    property_labels = []
    fmt = []
    regex = []
    if anc_label:
        property_labels.append("anc_label")
        if anc_link_path:
            fmt.append(anc_link_path)
        else:            
            fmt.append("{}")

    if event_label:
        property_labels.append("event")
        fmt.append("{}")

    if bs_label:
        property_labels.append("bs_support")
        fmt.append("{:.0f}")

    if pp_label:
        property_labels.append("anc_pp")
        fmt.append("{:.2f}")

    # Draw propertly labels on the tree
    if len(property_labels) > 0:
        pt.draw_node_labels(property_labels=property_labels,
                            fmt_string="|".join(fmt),
                            position=label_position,
                            position_offset=label_position_offset,
                            text_style=label_text_style)


    # Plot legends
    pt.draw_scale_bar()
    pt.draw_node_legend(label_renamer={"event":"events",
                                       "bs_support":"branch support",
                                       "anc_pp":"anc post. prob."})

    # Figure out plotting to file.
    if output_file is not None:
        pt.render(output_file)

    # If this is in a notebook, render it and return so it appears
    ret = None
    if topiary._in_notebook or return_canvas:
        ret = pt.canvas

    return ret
