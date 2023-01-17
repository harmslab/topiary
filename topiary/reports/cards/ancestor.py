"""
Functions for generating bootstrap cards with various types of information
on them.
"""

from topiary._private.supervisor import Supervisor
from topiary.draw.ancestor_data import plot_ancestor_data

from topiary.reports.elements import create_element
from topiary.reports.elements import create_card
from topiary.reports.elements import create_icon_row
from topiary.reports.elements import df_to_table
from topiary.reports.elements import sequence_box
from topiary.reports.elements import create_info_modal

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

import shutil
import os
import glob

overview_help_text = \
"""
The table holds general information about the ancestor. 

<b>Ancestor type</b> indicates the evolutionary event corresponding to the 
ancestral node:  "speciation", "duplication", or "transfer". It is only meaningful for 
reconciled trees. <b>Number of extant descendants</b> indicates the total number 
of modern sequences that arose from this ancestor. <b>Taxonomic distribution of
descendants</b> is the minimum taxonomic classification that encompasses the 
species of all descendants of this ancestor. An ancestor with only 
human and chimpanzee descendants would be "Hominini". An ancestor with human,
chimpanzee, and tarsier descendants would be "Simiiformes". It is only 
meaningful for reconciled trees. <b>Descendant paralog calls</b> is the 
the fraction of the descendants that were called as a particular paralog type
by reciprocal BLAST. <b>Mean posterior probability</b> is the mean of the 
posterior probability for the best reconstruction over all non-gap sites. It 
ranges from 1.0 (strong support at all sites) to 0.05 (completely ambiguous at
all sites). <b>Number of ambiguous sites</b> is the number of sites where the
posterior probability of the next-best reconstructed site is greater than alt_cutoff. 
<b>Number of ambiguous gaps</b> is the number of sites where it is unclear (by
maximum parsimony) if the position should be be reconstructed as an amino acid
or as a gap. <b>Branch support</b> is the branch support for the reconstructed
node. It ranges from 0 (no support) to 100 (high support).

The linked csv file has site-by-site statistics on the reconstruction. The 
<b>site_type column</b> will be one of: "good": the alt amino acid has posterior
probability (pp) below alt_cutoff; "ambiguous_similar": alt amino acid has pp
above alt_cutoff, but ML and alt amino acids are similar (e.g., T vs. S); 
"ambiguous_dissimilar": alt amino acid has pp above alt_cutoff, with ML and alt
dissimilar (e.g., E vs. F); "possible_gap": ambiguous whether this should be gap
or not; or "gap": is a gap. The <b>entropy column</b> is the Shannon entropy of the
posterior probabilities for all amino acids at that site. It ranges from 0 (one
amino acid has pp of 1; all others have pp of 0) to 3 (all twenty amino acids
have pp = 0.05). 

The linked fasta file holds the ML and altAll sequences with gaps removed.

The linked pdf file has the posterior probability plot. (Note: this pdf file is
likely easier to edit in Illustrator or Inkscape than the svg plot shown in the
report). 
"""

mlseq_help_text = \
"""
Maximum likelihood sequence of this ancestor. This sequence uses the amino acid
with the highest posterior probability at each position. The sequence is colored
by posterior probability. The gradient goes from <span style="color:red"><b>Red</b></span>,
indicating a posterior probability of 0.5 or less, to <b>Black</b>, indicating a 
posterior probability of 1.0.
"""

altall_help_text = \
"""
The altAll sequence starts with the maximum likelihood sequence, and then 
replaces ambiguous sites with the next-most-probable reconstructed amino acids. 
Ambiguous sites are defined as those where the posterior probability of the 
next-most-probable site is above alt_cutoff. (In the posterior probability plot,
these are sites where the red point is above the dashed line). One can determine
if an experimental result is robust to phylogenetic uncertainty by characterizing
both the maximum likelihood and altAll versions of an ancestor. In the graphic,
the maximum likelihood sites are white; the sites replaced with an alternate 
reconstruction are shown in black. For more information, see
<a href="https://doi.org/10.1093/molbev/msw223">Eick et al.</a>.
"""

pp_plot_help_text = \
"""
The <b>black</b> points show the posterior probability of the most likely amino
acid at each position. The distribution of these probabilities is given by the
<span style="color:gray;">gray</span> histogram on the right. The average posterior
probability is the mean of these values. The <span style="color:red;">red</span>
points show the probability of the second most likely amino acid at each position,
with its distribution in <span style="color:red;">red</span> on the right. The
horizontal dashed line indicates the alt_cutoff. Any
<span style="color:red;">red</span> point above this dashed line will be
included in the altAll sequence. Shaded <span style="color:gray;">gray</span>
regions indicate gaps; vertical <span style="color:purple;">purple</span> dashed
lines represent ambiguously gapped positions. The locations of these gaps can 
be identified in the csv file linked to in the ancestor overview. 
"""

def create_ancestor_card(anc_dict,
                         output_directory,
                         ancestor_directory,
                         event_color,
                         event_name):
    """
    Create a card holding all of the ancestors as individual accordion entries. 
    """

    anc_csv = glob.glob(os.path.join(ancestor_directory,
                                     "output",
                                     "*ancestors*",
                                     "ancestor-data.csv"))

    anc_csv = anc_csv[0]
    anc_dir = os.path.dirname(anc_csv)
    anc_df = pd.read_csv(anc_csv)
    alt_cutoff = Supervisor(ancestor_directory).run_parameters["alt_cutoff"]

    shutil.copy(anc_csv,output_directory)
    shutil.copy(os.path.join(anc_dir,"ancestors.fasta"),output_directory)

    anc_out = []
    start, _ = create_element("div",{"class":"accordion",
                                     "id":"ancAccordion"})
    anc_out.append(start)
    
    # Get sorted list of ancestors
    anc_seen = list(anc_dict.keys())
    anc_seen = [(int(a[3:]),a) for a in anc_seen]
    anc_seen.sort()
    anc_seen = [a[1] for a in anc_seen]
    for a in anc_seen:

        anc_id = f"a{a[3:]}"
    
        df = anc_df.loc[anc_df.anc == a,:]
        df.to_csv(os.path.join(output_directory,f"{a}.csv"),index=False)

        # Write out ML and altAll fasta file
        ml_seq = [s for s in df.ml_state if s != "-"]
        alt_seq = np.array(df.ml_state)
        alt_mask = np.array(df.alt_pp > alt_cutoff)
        alt_seq[alt_mask] = np.array(df.alt_state)[alt_mask]
        alt_seq = [s for s in alt_seq if s != "-"]

        f = open(os.path.join(output_directory,f"{a}.fasta"),"w")
        f.write(f">{a}|ML sequence\n")
        f.write("".join(ml_seq))
        f.write("\n")
        f.write(f">{a}|altAll sequence\n")
        f.write("".join(alt_seq))
        f.write("\n")
        f.close()

        fig, ax = plot_ancestor_data(df,
                                     width_ratio=6,
                                     alt_anc_pp=alt_cutoff,
                                     close_plot=False)
        fig.savefig(os.path.join(output_directory,f"{a}.svg"),bbox_inches = "tight")
        fig.savefig(os.path.join(output_directory,f"{a}_pp.pdf"),bbox_inches = "tight")

        plt.close(fig)

        taxonomic = anc_dict[a]["taxonomic_dist"]
        paralog_call = anc_dict[a]["paralog_call"]
        

        # ----------------------------------------------------------------------
        # Construct accordion item
        
        anc_out.append("<div class=\"accordion-item\">")

        # <header>
        anc_out.append("<div class=\"accordion-header\">")
        start, end = create_element("card",attributes={"class":"accordion-button",
                                                       "type":"button",
                                                       "data-bs-toggle":"collapse",
                                                       "data-bs-target":f"#collapse{anc_id}",
                                                       "aria-expanded":"true",
                                                       "aria-controls":f"collapse{anc_id}"})
        anc_out.append(start)
        if taxonomic is not None:
            anc_out.append(f"<h5 id=\"{anc_id}\">{a}: {taxonomic} {paralog_call}</h5>")
        else:
            anc_out.append(f"<h5 id=\"{anc_id}\">{a}: {paralog_call}</h5>")
        anc_out.append(end)
        anc_out.append("</div>")

        # </header>
        start, _   = create_element("div",attributes={"id":f"collapse{anc_id}",
                                                      "class":["accordion-collapse","collapse"],
                                                      "aria-labelledby":f"heading{anc_id}",
                                                      "data-bs-parent":"#ancAccordion"})

        anc_out.append(start)
        anc_out.append("<div class=\"accordion-body\">")
        
        # ----------------------------------------------------------------------
        # Create card holding general information for the ancestor
    
        event = anc_dict[a]["event"]
        e_color = event_color[event]
        e_name = event_name[event]
        event_html = f"<span style=\"color:{e_color}\">{e_name}</span>"

        paralogs = anc_dict[a]["paralogs"]

        mean_pp = float(anc_dict[a]["anc_pp"])
        bs_support = anc_dict[a]["bs_support"]
        if bs_support is None:
            bs_support = "N/A"
        else:
            bs_support = f"{int(round(float(bs_support),0))}"

        num_ambig_seq = np.sum(np.array(df.alt_pp > alt_cutoff))
        num_ambig_gap = np.sum(np.array(df.site_type == "possible gap"))
        
        stats_df = pd.DataFrame({"descriptions":["Ancestor type",
                                                 "Number of extant descendants",
                                                 "Taxonomic distribution of descendants",
                                                 "Descendant paralog calls",
                                                 "Mean posterior probability",
                                                 "Number of ambiguous sites",
                                                 "Number of ambiguous gaps",
                                                 "Branch support"],
                                "values":       [event_html,
                                                 anc_dict[a]["num_descendants"],
                                                 taxonomic,
                                                 paralogs,
                                                 f"{mean_pp:.2f}",
                                                 f"{num_ambig_seq}",
                                                 f"{num_ambig_gap}",
                                                 bs_support]})
        stats_html = df_to_table(stats_df,add_header=False,show_row_numbers=False)

        icon_html = create_icon_row([f"{a}.csv",f"{a}.fasta",f"{a}_pp.pdf"],
                                    [f"{a} csv",f"{a} fasta",f"{a} posterior probability plot"])

        help_html = create_info_modal(modal_text=overview_help_text,
                                      modal_title="Ancestor overview",
                                      extra_button_class="text-end")
    
        card_contents = "".join([stats_html,icon_html,help_html])
        
        stats_card = create_card(card_title=f"{a} overview",
                                 card_contents=card_contents,
                                 title_tag="h5")
        
        anc_out.append(stats_card)
        anc_out.append("<br/>")

        # ----------------------------------------------------------------------
        # Create card holding the ML sequence of the ancestor

        txt = np.array(df.ml_state)
        pp = np.array(df.ml_pp)
        ml_out = ["<p>Amino acids are colored by posterior probability: 0.5 (red); 1.0 (black)</p>"]
        ml_out.append(sequence_box(txt,
                                   prop_value=pp,
                                   prop_span=(0.5,1),
                                   color=["red","black"]))
        help_html = create_info_modal(modal_text=mlseq_help_text,
                                      modal_title="Maximum likelihood sequence",
                                      extra_button_class="text-end")
        help_html = f"<br/>{help_html}"
        ml_out.append(help_html)
        ml_card = create_card(card_title="Maximum likelihood sequence",
                              card_contents="".join(ml_out),
                              title_tag="h5")
        anc_out.append(ml_card)
        anc_out.append("<br/>")
        

        # ----------------------------------------------------------------------
        # Create card holding the altAll sequence of the ancestor
        
        txt = np.array(df.ml_state)
        mask = np.array(df.alt_pp > alt_cutoff)
        txt[mask] = np.array(df.alt_state)[mask]
        alt_out = ["<p>Black amino acids differ from the ML ancestor</p>"]
        alt_out.append(sequence_box(txt,
                                    prop_value=mask,
                                    prop_span=(0,1),
                                    color=["white","black"]))
     
        help_html = create_info_modal(modal_text=altall_help_text,
                                      modal_title="altAll sequence",
                                      extra_button_class="text-end")
        help_html = f"<br/>{help_html}"
        alt_out.append(help_html)

        alt_card = create_card(card_title="altAll sequence",
                               card_contents="".join(alt_out),
                               title_tag="h5")
        anc_out.append(alt_card)
        
        anc_out.append("<br/>")

        # ----------------------------------------------------------------------
        # Create card holding the posterior probability plot

        pp_plot_out = []
        pp_plot_out.append("<div class=\"text-center\">")
        pp_plot_out.append(f"<img src=\"{a}.svg\" alt=\"{anc_id} posterior probability plot\"/>")
        pp_plot_out.append("</div>")

        help_html = create_info_modal(modal_text=pp_plot_help_text,
                                      modal_title="posterior probability plot",
                                      extra_button_class="text-end")
        help_html = f"<br/>{help_html}"
        pp_plot_out.append(help_html)    


        pp_plot_card = create_card(card_title="Posterior probability plot",
                                   card_contents="".join(pp_plot_out),
                                   title_tag="h5")
        
        anc_out.append(pp_plot_card)
        anc_out.append("<br/>")                                              

        # ----------------------------------------------------------------------
        # Close out accordion item
        anc_out.append("</div>") # accordion-body
        anc_out.append("</div>") # accordion-collapse
        anc_out.append("</div>") # accordion-item

    # Accordion end
    anc_out.append("</div>") #ancAccordion div
    anc_html = "".join(anc_out)
    
    return create_card(card_contents=anc_html,
                       card_title="Ancestors",
                       title_tag="h4")

