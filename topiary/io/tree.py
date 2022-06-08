__description__ = \
"""
Read and write trees.
"""
__author__ = "Michael J. Harms"
__date__ = "2021-04-08"


import ete3
from ete3 import Tree
import dendropy as dp

def read_tree(tree,fmt=None):
    """
    Load a tree into an ete3 tree data structure.

    Parameters
    ----------
        tree: some sort of tree. can be an ete3.Tree (returns self), a dendropy
              Tree (converts to newick and drops root), a newick file or a
              newick string.
        fmt: format for reading tree from newick.  0-9 or 100. See ete3
             documentation for how these are read
             (http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees).
             As of ETE3.1.1, these numbers mean:


             |        ======  ==============================================
             |        FORMAT  DESCRIPTION
             |        ======  ==============================================
             |        0        flexible with support values
             |        1        flexible with internal node names
             |        2        all branches + leaf names + internal supports
             |        3        all branches + all names
             |        4        leaf branches + leaf names
             |        5        internal and leaf branches + leaf names
             |        6        internal branches + leaf names
             |        7        leaf branches + all names
             |        8        all names
             |        9        leaf names
             |        100      topology only
             |        ======  ==============================================

             if fmt is None, try to parse without a format descriptor, then these
             formats in numerical order.

    Return
    ------
        an ete3 tree object.
    """

    # Already an ete3 tree.
    if type(tree) is ete3.TreeNode:
        return tree

    # Convert dendropy tree into newick (drop root)
    if type(tree) is dp.Tree:
        tree = tree.as_string(schema="newick",suppress_rooting=True)

    # If we get here, we need to convert. If fmt is not specified, try to parse
    # without a format string.
    if fmt is None:


        try:
            t = Tree(tree)
        except ete3.parser.newick.NewickError:

            # Try all possible formats now, in succession
            w = "\n\nCould not parse tree without format string. Going to try different\n"
            w += "formats. Please check output carefully.\n\n"
            print(w)

            formats = list(range(10))
            formats.append(100)

            t = None
            for f in formats:
                try:
                    t = Tree(tree,format=f)
                    w = f"\n\nSuccessfully parsed tree with format style {f}.\n"
                    w += "Please see ete3 documentation for details:\n\n"
                    w += "http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees\n\n"
                    print(w)
                    break

                except ete3.parser.newick.NewickError:
                    continue

            if t is None:
                err = "\n\nCould not parse tree!\n\n"
                raise ValueError(err)

    else:
        # Try a conversion with the specified format
        t = Tree(tree,format=fmt)

    return t
