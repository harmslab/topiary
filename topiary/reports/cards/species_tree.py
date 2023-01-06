import topiary

from topiary.io.tree import read_tree

from topiary.reports.elements import canvas_to_html
from topiary.reports.elements import create_info_modal
from topiary.reports.elements import create_card
from topiary.reports.elements import create_icon_row

import shutil
import os

species_tree_help_text = \
"""
Species tree used for the gene/species tree reconciliation. This is a cladogram
(the branch lengths are not meaningful). The species names are aligned on the 
right, with dashed lines connecting the labels to the relevant branches.  This
tree is the latest synthetic tree downloaded from the 
<a href="https://tree.opentreeoflife.org/">Open Tree of Life</a> database. The
tree can be downloaded in newick format or as a pdf using the links to the left.
"""

def create_species_tree_card(supervisor,output_directory):

    species_T = read_tree(supervisor.species_tree)

    ott_to_name = dict(zip(supervisor.df.ott,supervisor.df.species))
    for n in species_T.traverse():
        if n.is_leaf():
            n.species = ott_to_name[n.name]

    shutil.copy(supervisor.species_tree,
                os.path.join(output_directory,"species-tree.newick"))

    tree_canvas = topiary.draw.species_tree(species_T,
                                            output_file=os.path.join(output_directory,"species-tree.pdf"),
                                            return_canvas=True)
    
    tree_html = canvas_to_html(tree_canvas)

    icon_html = create_icon_row(["species-tree.newick","species-tree.pdf"],
                                ["species tree as newick","species tree as a pdf"])

    help_html = create_info_modal(modal_text=species_tree_help_text,
                                  modal_title="Species tree",
                                  extra_button_class="text-end")
    help_html = f"<br/>{help_html}"
    tree_html = f"{tree_html}{icon_html}{help_html}"

    return create_card(card_title="Species tree for gene/species tree reconciliation",
                       card_contents=tree_html,title_tag="h4")


