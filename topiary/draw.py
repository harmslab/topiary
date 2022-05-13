__description__ = \
"""
Functions for plotting ete3 trees.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-05-12"

import topiary
import ete3

def species_tree(species_tree,output=None):
    """
    species_tree: ete3 Tree object holding a species tree
    output: None or output file. If None and work being done
            in a notebook, render inline. If not in a notebook,
            write to species_tree.pdf.
    """

    T = species_tree.copy()

    ts = ete3.TreeStyle()
    ts.draw_guiding_lines = True
    ts.guiding_lines_type = 0
    ts.show_scale = False
    ts.scale = 20

    for n in T.traverse():

        ns = ete3.NodeStyle()
        ns["size"] = 0
        ns["hz_line_width"] = 2
        ns["vt_line_width"] = 2

        n.set_style(ns)

        if n.is_leaf():
            n.name = ""

            txt = ete3.TextFace(n.species[0],fsize=14)
            txt.margin_left = 4
            n.add_face(txt,0,position="branch-right")


    if output is None:
        if topiary._in_notebook:
            return T.render("%%inline",tree_style=ts)

        else:
            output = "species_tree.pdf"

    T.render(output,tree_style=ts)
