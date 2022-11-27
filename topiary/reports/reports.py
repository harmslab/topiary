"""
Generate an html report for a topiary tree/ancestor inference.
"""

import topiary

from topiary._private import Supervisor

from topiary.draw import plot_ancestor_data

from topiary.reports.elements import create_output_directory
from topiary.reports.elements import canvas_to_html
from topiary.reports.elements import create_main_html
from topiary.reports.elements import df_to_table
from topiary.reports.elements import create_card
from topiary.reports.elements import create_element
from topiary.reports.elements import create_icon_row

from topiary.reports.quality import check_duplication

import pandas as pd
import numpy as np

import os
import glob
import shutil

EVENT_COLOR = {"D":"#64007F","L":"#BAD316","T":"#407E98","S":"#023E55"}
EVENT_NAME = {"S":"speciation","D":"duplication","L":"loss","T":"transfer"}

def _find_directory(calculation_directory):
    """
    Figure out the directory to use for the calculation. 

    Parameters
    ----------
    calculation_directory : str
        calculation directory specified by the user. Should either be a pipeline
        directory containing a bunch of individual calculations OR a single
        calculation directory (with run_parameters.json). If a pipeline 
        directory is specified, return the latest calculation and the ancestor
        directory if found. 
    
    Returns
    -------
    calculation_directory : str or None
        return path to calculation directory if found, None otherwise
    ancestor_directory : str or None
        return path to ancestor directory if found, None otherwise

    """

    if os.path.isfile(os.path.join(calculation_directory,"run_parameters.json")):
        return calculation_directory, None

    calc_dirs = []
    for c in os.listdir(calculation_directory):
        this_dir = os.path.join(calculation_directory,c)
        if os.path.isdir(this_dir):
            if os.path.exists(os.path.join(this_dir,"run_parameters.json")):
                sv = Supervisor(this_dir)
                if "ancestors" in sv.calc_type:
                    is_anc_dir = True
                else:
                    is_anc_dir = False

                try:
                    t = sv.run_parameters["creation_time"]
                except KeyError:
                    t = -1

                calc_dirs.append((t,this_dir,is_anc_dir))
    
    # Sort from oldest to newest
    calc_dirs.sort()
    calc_dirs = calc_dirs[::-1]

    # If there is a calculation in this directory
    if len(calc_dirs) > 0:

        # See if we can find an ancestor directory
        ancestor_dir = None
        for c in calc_dirs:

            # Is an ancestor directory
            if c[2]:
                ancestor_dir = c[1]
                break
        
        # Return calculation directory and ancestor directory if found
        return calc_dirs[0][1], ancestor_dir

    # No calculation found, no ancestor directory found
    return None, None

def _create_ancestor_card(anc_dict,
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
                                    alt_anc_pp=alt_cutoff)
        fig.savefig(os.path.join(output_directory,f"{a}.svg"),bbox_inches = "tight")
        fig.savefig(os.path.join(output_directory,f"{a}_pp.pdf"),bbox_inches = "tight")

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
        anc_out.append(f"<h5 id=\"{anc_id}\">{a}: {taxonomic} {paralog_call}</h5>")
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
        
        stats_df = pd.DataFrame({"descriptions":["Ancestor type",
                                                "Number of extant descendants",
                                                "Taxonomic distribution of descendants",
                                                "Descendant paralog calls",
                                                "Mean posterior probability",
                                                "Branch support"],
                                "values":      [event_html,
                                                anc_dict[a]["num_descendants"],
                                                taxonomic,
                                                paralogs,
                                                f"{mean_pp:.2f}",
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
        ml_out.append(topiary.reports.elements.sequence_box(txt,
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
        alt_out.append(topiary.reports.elements.sequence_box(txt,
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


def create_report(calculation_directory,
                  output_directory,
                  ancestor_directory=None,
                  html_description="Topiary report",
                  html_title="Topiary report",
                  overwrite=True):
    """
    Create a sharable html file and directory holding the results of a topiary
    calculation.

    Parameters
    ----------
    calculation_directory : str
        path to a directory with a calculation. This can either point to a 
        pipeline directory (e.g., output of ali_to_anc) OR a single calculation
        (e.g., output of generate_ml_tree). If a pipeline directory is 
        specified, the function will generate a report covering all aspects of
        the pipeline. 
    output_directory : str
        directory in which to write the report
    ancestor_directory : str, optional
        directory holding ancestors. If specified, this will override ancestors
        found automatically in calculation_directory. 
    overwrite : bool, default = False
        whether or not to overwrite an existing report. 
    """


    # Get the calculation directory
    calculation_directory, inferred_anc_dir = _find_directory(calculation_directory)
    if ancestor_directory is None:
        ancestor_directory = inferred_anc_dir

    # Load calculation into supervisor
    supervisor = Supervisor(calculation_directory)
    
    # Get information about ancestors (pp support, evolutionary event, 
    # bs_support, etc.) Stick into anc_dict
    T = topiary.draw.core.load_trees(supervisor.output_dir,
                                     prefix=supervisor.tree_prefix)

    if "recip_paralog" in supervisor.df.columns:
        p_column = "recip_paralog"
    else:
        p_column = "name"

    anc_dict = {}
    for n in T.traverse():

        if not n.is_leaf():

            if n.anc_label is not None:

                anc_name = f"anc{n.anc_label[1:]}"
                anc_dict[anc_name] = {}
                
                for f in n.features:
                    try:
                        anc_dict[anc_name][f] = n.__dict__[f]
                    except KeyError:
                        anc_dict[anc_name][f] = n.__dict__[f"_{f}"]

                # Get leaves for this ancestor
                leaves = [l.name for l in n.get_leaves()]
                anc_dict[anc_name]["num_descendants"] = len(leaves)

                # Get dataframe holding descendants of this ancestor
                mask = supervisor.df.uid.isin(leaves)
                descendant_df = supervisor.df.loc[mask,:]

                # Get name of taxonomic grouping for the ancestors
                mrca = topiary.opentree.ott_to_mrca(ott_list=descendant_df.ott)
                anc_dict[anc_name]["taxonomic_dist"] = mrca["ott_name"]

                # Get paralog descendant percentages
                paralogs = np.unique(descendant_df.loc[:,p_column])
                paralogs = [(np.sum(descendant_df.loc[:,p_column] == p),p) for p in paralogs]
                paralogs = [p for p in paralogs if p[0] > 0]
                paralogs.sort()
                paralogs = paralogs[::-1]
                total = sum([p[0] for p in paralogs])
                paralogs = [(p[0]/total,p[1]) for p in paralogs]

                p_str = [f"{p[1]}: {100*p[0]:.1f}%" for p in paralogs]
                anc_dict[anc_name]["paralogs"] = ",".join(p_str)

                anc_dict[anc_name]["paralog_call"] = "|".join([p[1] for p in paralogs if p[0] > 0.1])
                
    
    # Build output directory
    create_output_directory(output_directory=output_directory,
                            overwrite=overwrite)

    supervisor.df.to_csv(os.path.join(output_directory,"dataframe.csv"))

    # Card stack will hold all of the generated html output
    card_stack = []

    # -------------------------------------------------------------------------
    # Title card

    title_html = create_card("Topiary output",calculation_directory,"h3")
    card_stack.append(title_html)

    # -------------------------------------------------------------------------
    # warning card

    if supervisor.tree_prefix == "reconciled":
        expect, obs, duplication_df = check_duplication(supervisor,T,p_column)
        if obs > expect > 0:

            out = []
            out.append(f"We expected {expect} duplications,")
            out.append(f"but observed {obs} duplications.")
            out.append("These extra duplications could be real, but could")
            out.append("also reflect model violation or incomplete lineage")
            out.append("sorting. The unexpected duplications--and number")
            out.append("of affected descendants--are shown below.")

            txt = " ".join(out)
            out = f"<p>{txt}</p><br/>"
            warning_table = df_to_table(duplication_df,show_row_numbers=False)
            warning_txt = f"{out}{warning_table}"

            warning_html = create_card("<span color=\"var(--red-2)\">Warning</span>",warning_txt,"h5")
            card_stack.append(warning_html)


    # -------------------------------------------------------------------------
    # Input card

    this_df = supervisor.df.loc[supervisor.df.keep,:]
    num_input = len(this_df)
    combined_paralogs = list(np.unique(this_df.loc[:,p_column]))
    combined_paralogs = [p.split("|") for p in combined_paralogs]
    paralogs = []
    for c in combined_paralogs:
        paralogs.extend(c)
    paralogs = np.array(np.unique(paralogs))
    paralogs.sort()
    num_paralogs = len(paralogs)
    mrca = topiary.opentree.ott_to_mrca(ott_list=np.unique(this_df.ott))
    taxonomic_distribution = mrca["ott_name"]

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

    card_stack.append(input_html)
    
    # -------------------------------------------------------------------------
    # Run parameters

    model = supervisor.model
    if supervisor.tree_prefix == "reconciled":
        reconciled = True
    else:
        reconciled = False

    ## FIGURE OUT IF WE GOT SUPPORTS SOMEWHERE XX
    supports = "XX"

    param_df = pd.DataFrame({"name":["Evolutionary model",
                                     "Reconciled gene and species tree",
                                     "Supports calculated"],
                             "value":[model,reconciled,supports]})
    param_table = df_to_table(param_df,add_header=False,show_row_numbers=False)
    param_html = create_card("Run parameters",card_contents=param_table)
    card_stack.append(param_html)

    # -------------------------------------------------------------------------
    # Icon card

    icon_html = create_icon_row([f"ancestor-data.csv",f"ancestors.fasta",f"summary-tree.pdf"],
                                [f"all ancestors csv",f"all ancestors fasta",f"ancestor tree pdf"])

    card_stack.append(create_card(card_title="Summary files",
                                  card_contents=icon_html))

    # -------------------------------------------------------------------------
    # Tree card

    if ancestor_directory is not None:
        anc_link_path = "<a href=\"#{anc_label}\">{anc_label}</a>"
    else:
        anc_link_path = None

    tree_canvas = topiary.draw.tree(supervisor,
                                    anc_link_path=anc_link_path,
                                    return_canvas=True)
    
    tree_html = canvas_to_html(tree_canvas)

    card_stack.append(create_card(card_contents=tree_html))

    # -------------------------------------------------------------------------
    # Ancestors card

    if ancestor_directory is not None:
        anc_card = _create_ancestor_card(anc_dict,
                                         output_directory,
                                         ancestor_directory,
                                         event_color=EVENT_COLOR,
                                         event_name=EVENT_NAME)
        card_stack.append(anc_card)
        


    # Assemble final report html

    container = ["<div class=\"container-lg px-4 gy-4\">"]
    container.append("".join(card_stack))
    container.append("</div>")
    container_html = "".join(container)

    top, bottom = create_main_html(description=html_description,
                                   title=html_title)

    out = [top,container_html,bottom]

    f = open(os.path.join(output_directory,"index.html"),"w")
    f.write("\n".join(out))
    f.close()

    


