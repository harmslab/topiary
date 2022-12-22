"""
Functions for generating bootstrap cards with various types of information
on them.
"""

import topiary

from topiary._private.supervisor import Supervisor
from topiary.draw.ancestor_data import plot_ancestor_data
from topiary.reports.elements import create_element
from topiary.reports.elements import create_card
from topiary.reports.elements import create_icon_row
from topiary.reports.elements import df_to_table
from topiary.reports.elements import sequence_box

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

import shutil
import os
import glob

def create_ancestor_card(anc_dict,
                         output_directory,
                         ancestor_directory,
                         event_color,
                         event_name):
    """
    Create a bootstrap card holding all of the ancestors as individual 
    accordion entries. 
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
    shutil.copy(os.path.join(anc_dir,"..","summary-tree.pdf"),output_directory)

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
        df.to_csv(os.path.join(output_directory,f"{a}.csv"))

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

        mean_pp = np.float(anc_dict[a]["anc_pp"])
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

        card_contents = "".join([stats_html,icon_html])
        
        stats_card = create_card(card_title=f"{a} overview",
                                 card_contents=card_contents)
        
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
        ml_card = create_card(card_title="Maximum likelihood sequence",
                            card_contents="".join(ml_out))
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
     
        alt_card = create_card(card_title="altAll sequence",
                            card_contents="".join(alt_out))
        anc_out.append(alt_card)
        anc_out.append("<br/>")

        # ----------------------------------------------------------------------
        # Create card holding the posterior probability plot

        pp_plot_out = []
        pp_plot_out.append("<div class=\"text-center\">")
        pp_plot_out.append(f"<img src=\"{a}.svg\" alt=\"{anc_id} posterior probability plot\"/>")
        pp_plot_out.append("</div>")
    
        pp_plot_card = create_card(card_title="Posterior probability plot",
                                card_contents="".join(pp_plot_out))
        
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
    
    return create_card(card_contents=anc_html)


def create_input_card(supervisor,p_column=None):
    """
    Create a card describing the input to the analysis. 
    """
    
    # -------------------------------------------------------------------------
    # Input card

    this_df = supervisor.df.loc[supervisor.df.keep,:]

    if p_column is None:
        if "recip_paralog" in this_df.columns:
            p_column = "recip_paralog"
        else:
            p_column = "name"

    num_input = len(this_df)
    combined_paralogs = list(np.unique(this_df.loc[:,p_column]))
    combined_paralogs = [p.split("|") for p in combined_paralogs]
    paralogs = []
    for c in combined_paralogs:
        paralogs.extend(c)
    paralogs = np.array(np.unique(paralogs))
    paralogs.sort()
    num_paralogs = len(paralogs)
    
    if np.sum(pd.isnull(this_df.ott)) == 0:
        ott_list = np.unique(this_df.ott)
        mrca = topiary.opentree.ott_to_mrca(ott_list=ott_list)
        taxonomic_distribution = mrca["ott_name"]
    else:
        taxonomic_distribution = None

    input_df = pd.DataFrame({"name":["Number of sequences",
                                     "Number of paralogs",
                                     "Paralogs",
                                     "Taxonomic distribution"],
                           "values":[num_input,
                                     num_paralogs,
                                     ",".join(paralogs),
                                     taxonomic_distribution]})

    input_table = df_to_table(input_df,add_header=False,show_row_numbers=False)                                 
    input_icons = create_icon_row(files_to_link=["dataframe.csv"],
                                  descriptions=["input dataframe"])
    input_card_input = f"<div>{input_table}{input_icons}</div>"

    input_html = create_card("Input",card_contents=input_card_input)

    return input_html


def create_model_card(supervisor,
                      output_directory):

    best_model = supervisor.model

    model_directory = supervisor.calc_dir
    model_comp = pd.read_csv(os.path.join(model_directory,"output","model-comparison.csv"))
    model_comp = model_comp.loc[:, ~model_comp.columns.str.contains('^Unnamed')]
    num_combos = len(model_comp)

    # Copy model comparison into the output directory
    shutil.copy(os.path.join(model_directory,"output","model-comparison.csv"),
                os.path.join(output_directory,"model-comparison.csv"))

    # make html for the top ten models
    top_ten_models = df_to_table(model_comp.iloc[:10,:],show_row_numbers=False,float_fmt="{:.3f}")

    s3, e3 = create_element("h5")
    s4, e4 = create_element("h5")

    out = []
    out.append(f"{s3}Best model: {best_model}{e3}")
    out.append(f"{s4}Number of models: {num_combos}{e4}")


    con_s, con_e = create_element("div",{"class":["text-center"]})
    rs, re = create_element("div",{"class":["row"]})
    cs, ce = create_element("div",{"class":["col"]})

    out.append(f"{con_s}{rs}{cs}{top_ten_models}{ce}")

    key_dict = {"value":["model name",
                        "log likelihood",
                        "Akaike Information Criterion",
                        "Akaike Information Criterion (corrected for sampling)",
                        "Bayesian Information Criterion",
                        "Number of parameters",
                        "model posterior probability"],
                "key":["model","L","AIC","AICc","BIC","N","p"]}
                        
    key_html = df_to_table(pd.DataFrame(key_dict),add_header=False,show_row_numbers=False)
    icon_html = create_icon_row(["model-comparison.csv"],["Table comparing models"])
    #out.append(f"{s4}Table key{e4}")

    out.append(f"{cs}{key_html}<br/>{icon_html}{ce}")
    out.append(f"{re}{con_e}")

    

    #out.append(icon_html)

    html = "".join(out)

    html = create_card(card_title="Evolutionary model selection",
                       card_contents=html,
                       title_tag="h4")
    
    return html
