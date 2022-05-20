__description__ = \
"""
Convenience function for user. User can draw a tree for any valid calcualtion
type using a single function call.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-05-19"

import topiary
from topiary.external.interface import read_previous_run_dir
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
        run_dir: directory containing an ML tree run
        output_file: output file to write tree. If running in a notebook, will
                     return to notebook and write to this file. If not in a notebook
                     and not specified, will write out ml-tree.pdf.
        tip_columns: label the tree tips as "|".join(tip_columns). If tip_columns
                     is ["species","nickname"], tips will have names like
                     'Homo sapiens|LY96'.
        tip_name_separator: string to separate columns in tip names ("|" in
                            tip_columns example above.)
        fontsize: fontsize in points for tip labels and support values
        circle_radius: circle size for internal nodes (fraction of total width)
        value_range: tuple holding minimum and maximum values for support color map.
                     (0,100) means miminum support value is 0, maximum support value
                     is 100. [Only relevant for ML tree and ancestor
                     calculations, where it will map to support values or
                     posterior probabilities, respectively.]
        R: tuple holding red channel value for lowest support and highest
           support. [Only relevant for ML tree and ancestor calculations.]
        G: tuple holding green channel value for lowest support and highest
           support. [Only relevant for ML tree and ancestor calculations.]
        B: tuple holding blue channel value for lowest support and highest
           support. [Only relevant for ML tree and ancestor calculations.]
        event_color_map: map between events (duplication (D), loss (L),
                         speciation (S), or transfer (T) and node color).[Only
                         relevant for a reconcilation tree].
        width: width of total tree (pixels)
        space_between_taxa: number of pixels between taxa, sets hight.
        line_width: width of tree lines (pixels)
        df: dataframe (overides whatever is in run_dir/output/dataframe.csv)

    Return
    ------
        None or rendered output if in notebook
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
