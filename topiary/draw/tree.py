"""
Draw tree output from a maximum likelihood tree calculation, gene/species tree
reconcilation, or ancestral sequence reconstruction. 
"""

import topiary
from topiary.external._interface import read_previous_run_dir
import inspect

def tree(run_dir,
         output_file=None,
         tip_columns=["species","nickname"],
         tip_name_separator="|",
         fontsize=36,
         circle_radius=0.025,
         value_range=None,
         R=(0.95,0),G=(0.95,0),B=(1,1),
         event_color_map={"D":"green","L":"red","S":"blue","T":"pink"},
         width=800,space_between_taxa=10,line_width=4,
         df=None):
    """
    Draw tree output from an ML tree calculation, gene/species tree
    reconcilation, or ancestral sequence reconstruction. Plot type is determined
    based on the type of calculation passed in as run_dir.

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
    event_color_map : dict
        map between evolutionary events (D,L,S,T) and named colors.
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

    prev_run = read_previous_run_dir(run_dir)
    calc_type = prev_run["calc_type"]

    if calc_type == "ml_tree":

        # If no value_range specified, get the default from ml_tree
        if value_range is None:
            sig = inspect.signature(topiary.draw.ml_tree)
            value_range = sig.parameters["value_range"].default

        return topiary.draw.ml_tree(run_dir=run_dir,
                                    output_file=output_file,
                                    tip_columns=tip_columns,
                                    tip_name_separator=tip_name_separator,
                                    fontsize=fontsize,
                                    circle_radius=circle_radius,
                                    value_range=value_range,
                                    R=R,G=G,B=B,
                                    width=width,
                                    space_between_taxa=space_between_taxa,
                                    line_width=line_width,
                                    df=df)

    elif calc_type == "reconciliation":

        return topiary.draw.reconciliation_tree(run_dir=run_dir,
                                                output_file=output_file,
                                                tip_columns=tip_columns,
                                                tip_name_separator=tip_name_separator,
                                                fontsize=fontsize,
                                                circle_radius=circle_radius,
                                                event_color_map=event_color_map,
                                                width=width,
                                                space_between_taxa=space_between_taxa,
                                                line_width=line_width,
                                                df=df)

    elif calc_type == "ancestors":

        # If no value_range specified, get the default from ml_tree
        if value_range is None:
            sig = inspect.signature(topiary.draw.ancestor_tree)
            value_range = sig.parameters["value_range"].default

        return topiary.draw.ancestor_tree(run_dir=run_dir,
                                          output_file=output_file,
                                          tip_columns=tip_columns,
                                          tip_name_separator=tip_name_separator,
                                          fontsize=fontsize,
                                          circle_radius=circle_radius,
                                          value_range=value_range,
                                          R=R,G=G,B=B,
                                          width=width,
                                          space_between_taxa=space_between_taxa,
                                          line_width=line_width,
                                          df=df)

    else:
        err = f"\nCould not draw a tree for a calculation of type '{calc_type}'\n\n"
        raise ValueError(err)
