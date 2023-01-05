import topiary

from topiary.reports.elements import create_card
from topiary.reports.elements import canvas_to_html
from topiary.reports.elements import create_info_modal
from topiary.reports.elements import create_icon_row

from topiary.io.tree import write_trees
from topiary.draw.core import create_name_dict

import os

tree_plot_help_text = \
"""
Branch lengths represent the average number of amino acid substitutions per site
on that branch and can be estimated using the scale bar. The tree is rooted by
the midpoint method (for the gene tree) or using the species tree (reconciled
tree). 

Reconstructed ancestral sequences at each node are labeled with a unique name
(a1, a2, etc.). These labels are links to more detailed information about that
ancestor. 

Each ancestral node is annotated with up to three concentric circles. The
outermost circle denotes the evolutionary event for that node. (This is only 
relevant for reconciled trees). Duplication events are marked in 
<span style="color:#64007F">purple</span> and transfers in
<span style="color:#407E98">teal</span>. Speciation events are not labeled. 
The middle circle indicates bootstrap support, ranging from weak (50 or below;
white) to strong (100; <b>black</b>). The center circle indicates average ancestor
posterior probability across all non-gap sites, ranging from weak (0.7 or below;
white) to strong (1.0; <span style="color:#DC801A">orange</span>). 

The newick file holds up to four phylogenetic trees in a single file. The tips
are labeled with species and paralog. The internal nodes are labeled with
ancestor names (a1, a2, etc.), ancestor posterior probabilities (0.0-1.0), 
branch supports (0-100), and, for reconciled trees, evolutionary events 
(D: duplication, T: transfer, S: speciation, L: loss). This file should be 
readable by all phylogenetic tree plotting and editing software. 

The pdf holds the phylogenetic tree. In our experience, the pdf version of the
tree is better than the svg file shown in the main panel for loading into Adobe
Illustrator or Inkscape for figure creation. 
"""

def create_asr_tree_card(supervisor,output_directory,ancestor_directory,T):

    if ancestor_directory is not None:
        anc_link_path = "<a href=\"#{anc_label}\">{anc_label}</a>"
    else:
        anc_link_path = None

    # Figure out what to call tips of tree
    if "recip_paralog" in supervisor.df.columns:
        tip_columns = ["species","recip_paralog","uid"]
    elif "nickname" in supervisor.df.columns:
        tip_columns = ["species","nickname","uid"]
    else:
        tip_columns = ["species","name","uid"]

    # Create dictionary mapping between uid and pretty name format for tips
    name_dict = create_name_dict(df=supervisor.df,
                                 tip_columns=tip_columns,
                                 disambiguate=True)

    # Write out trees as a newick file
    write_trees(T,
                name_dict=name_dict,
                out_file=os.path.join(output_directory,"asr-trees.newick"))

    # -------------------------------------------------------------------------
    # Draw tree and write as pdf

    tree_canvas = topiary.draw.tree(supervisor,
                                    anc_link_path=anc_link_path,
                                    output_file=os.path.join(output_directory,
                                                             "asr-tree.pdf"),
                                    return_canvas=True)
    
    tree_html = canvas_to_html(tree_canvas)

    icon_html = create_icon_row(["asr-trees.newick","asr-tree.pdf"],
                                ["asr trees as a single newick",
                                 "asr tree as a pdf"])

    help_html = create_info_modal(modal_text=tree_plot_help_text,
                                  modal_title="Tree plot",
                                  extra_button_class="text-end")
    help_html = f"<br/>{help_html}"
    tree_html = f"{tree_html}{icon_html}{help_html}"

    return create_card(card_title="Tree used for ASR",card_contents=tree_html,title_tag="h4")