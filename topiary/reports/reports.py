"""
Generate various html reports for a topiary tree/ancestor inferences.
"""

import topiary

from topiary._private import Supervisor

from topiary.io.tree import load_trees


from topiary.reports.elements import create_output_directory
from topiary.reports.elements import canvas_to_html
from topiary.reports.elements import create_main_html
from topiary.reports.elements import df_to_table
from topiary.reports.elements import create_card
from topiary.reports.elements import create_row
from topiary.reports.elements import create_icon_row
from topiary.reports.elements import create_info_modal

from topiary.reports.cards import create_ancestor_card
from topiary.reports.cards import create_input_card
from topiary.reports.cards import create_model_card
from topiary.reports.cards import create_duplications_card
from topiary.reports.cards import create_species_tree_card
from topiary.reports.cards import create_asr_tree_card
from topiary.reports.cards import create_param_card

import pandas as pd
import numpy as np

import os
import shutil

EVENT_COLOR = {"D":"#64007F","L":"#BAD316","T":"#407E98","S":"#023E55",None:"#000000"}
EVENT_NAME = {"S":"speciation","D":"duplication","L":"loss","T":"transfer",None:'N/A'}

summary_files_help_text = \
"""
The linked csv file has site-by-site statistics on all reconstructed ancestors. The
<b>site_type column</b> will be one of: "good": the alt amino acid has posterior
probability (pp) below alt_cutoff; "ambiguous_similar": alt amino acid has pp
above alt_cutoff, but ML and alt amino acids are similar (e.g., T vs. S); 
"ambiguous_dissimilar": alt amino acid has pp above alt_cutoff, with ML and alt
dissimilar (e.g., E vs. F); "possible_gap": ambiguous whether this should be gap
or not; or "gap": is a gap. The <b>entropy column</b> is the Shannon entropy of the
posterior probabilities for all amino acids at that site. It ranges from 0 (one
amino acid has pp of 1; all others have pp of 0) to 3 (all twenty amino acids
have pp = 0.05). 

The fasta file holds all maximum likelihood and altAll ancestors.
"""

gene_tree_help_text = \
"""
Ancestors calculated using the maximum likelihood gene tree rather than the 
reconciled gene species tree.
"""

reconciled_tree_help_text = \
"""
Ancestors calculated using the reconciled gene/species tree rather than the 
maximum likelihood gene tree.
"""


def _find_directories(calculation_directory):
    """
    Figure out the directories in a pipeline to use for a report. 

    Parameters
    ----------
    calculation_directory : str
        calculation directory specified by the user. Should either be a pipeline
        directory containing a bunch of individual calculations OR a single
        calculation directory (with run_parameters.json). 
    
    Returns
    -------
    calc_dirs : dict
        dictionary with keys "model", "gene" and "reconciled" pointing to the
        results of a pipeline calculation. For "model", the value is the 
        directory for the model inference. For the gene or reconciled keys, 
        the values are dictionaries with keys "anc" and "tree" pointing to the
        last completed directory with an ancestor or tree inference, respectively.
    """

    calc_dirs = {"model"     :None,
                 "gene"      :{"anc":[],
                               "tree":[]},
                 "reconciled":{"anc":[],
                               "tree":[]}}

    if os.path.isfile(os.path.join(calculation_directory,"run_parameters.json")):
        input_calc_dirs = [calculation_directory]
    else:
        input_calc_dirs = [os.path.join(calculation_directory,c) for c in 
                           os.listdir(calculation_directory)]

    for c in input_calc_dirs:

        if os.path.isdir(c):
            if os.path.exists(os.path.join(c,"run_parameters.json")):
                sv = Supervisor(c)

                # Run has not finished
                if sv.status != "complete":
                    continue

                # creation time
                try:
                    t = sv.run_parameters["creation_time"]
                except KeyError:
                    t = -1                

                if sv.calc_type == "find_best_model":
                    calc_dirs["model"] = c

                # tree type (will be gene, reconciled, or None)
                tree_type = sv.tree_prefix
                if tree_type is None:
                    continue

                # Is an ancestor dir, append to both anc and tree
                if "ancestors" in sv.calc_type:
                    calc_dirs[tree_type]["anc"].append((t,c))
                    calc_dirs[tree_type]["tree"].append((t,c))

                # If a tree dir, append ot only anc
                else:
                    calc_dirs[tree_type]["tree"].append((t,c))

    for tree_type in ["gene","reconciled"]:

        # Take last anc and tree calcs seen
        for type_key in ["anc","tree"]:

            # no calc -- set to None
            if len(calc_dirs[tree_type][type_key]) == 0:
                calc_dirs[tree_type][type_key] = None
                continue

            calc_dirs[tree_type][type_key].sort()
            calc_dirs[tree_type][type_key] = calc_dirs[tree_type][type_key][-1][1]


    return calc_dirs


def tree_report(tree_directory,
                output_directory,
                ancestor_directory=None,
                html_description="Topiary report",
                html_title="Topiary report",
                overwrite=False,
                create_zip_file=True):
    """
    Create a sharable html file and directory holding the results of a topiary
    calculation.

    Parameters
    ----------
    tree_directory : str
        path to a directory with a tree calculation
    output_directory : str
        directory in which to write the report
    ancestor_directory : str, optional
        directory holding ancestors that match the tree
    overwrite : bool, default = False
        whether or not to overwrite an existing report. 
    create_zip_file : bool, default = True
        whether or not to create a zip file of the report directory. This zip
        file will be self-contained and is a way to share the output of the 
        calculation with someone else. 
    """

    # Load calculation into supervisor
    supervisor = Supervisor(tree_directory)
    
    # Load trees in as an ete3 tree
    T = load_trees(supervisor.output_dir,
                   prefix=supervisor.tree_prefix)

    # Figure out column with paralog calls
    if "recip_paralog" in supervisor.df.columns:
        p_column = "recip_paralog"
    else:
        p_column = "name"

    # Get information about ancestors (pp support, evolutionary event, 
    # bs_support, etc.) Stick into anc_dict
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
                if supervisor.tree_prefix == "reconciled":
                    mrca = topiary.opentree.ott_to_mrca(ott_list=descendant_df.ott)
                    anc_dict[anc_name]["taxonomic_dist"] = mrca["ott_name"]
                else:
                    anc_dict[anc_name]["taxonomic_dist"] = None

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


    print(f"Generating report in {output_directory}",flush=True)
    
    # -------------------------------------------------------------------------
    # Create output directory and put in dataframe

    # Build output directory
    create_output_directory(output_directory=output_directory,
                            overwrite=overwrite)

    supervisor.df.to_csv(os.path.join(output_directory,"dataframe.csv"),index=False)

    # Card stack will hold all of the generated html output
    card_stack = []

            
    # -------------------------------------------------------------------------
    # Title card

    title_html = create_card(html_description,tree_directory,"h3")
    card_stack.append(f"{title_html}<br/>")

    # -------------------------------------------------------------------------
    # Warning card

    if supervisor.tree_prefix == "reconciled":
        duplications_card = create_duplications_card(supervisor,T,p_column)
        card_stack.append(duplications_card)

    # -------------------------------------------------------------------------
    # input card

    input_html = create_input_card(supervisor,p_column)
    
    # -------------------------------------------------------------------------
    # Run parameters

    param_html = create_param_card(supervisor,anc_dict,ancestor_directory)

    # -------------------------------------------------------------------------
    # Icon card

    icon_html = create_icon_row(["ancestor-data.csv","ancestors.fasta"],
                                ["all ancestors csv","all ancestors fasta"])

    help_html = create_info_modal(modal_text=summary_files_help_text,
                                  modal_title="Summary files",
                                  extra_button_class="text-end")
    help_html = 6*"<br/>"+ help_html
    icon_html = f"<br/>{icon_html}{help_html}"

    icon_html = create_card(card_title="Summary files",
                            card_contents=icon_html,title_tag="h4")

    # Combine input, param, and icon into a single row
    card_stack.append(create_row([input_html,param_html,icon_html]))

    # -------------------------------------------------------------------------
    # Species tree card

    if supervisor.tree_prefix == "reconciled":
        species_tree_card = create_species_tree_card(supervisor,output_directory)
        card_stack.append(species_tree_card)

    # -------------------------------------------------------------------------
    # ASR tree card

    asr_tree_card = create_asr_tree_card(supervisor,output_directory,ancestor_directory,T)
    card_stack.append(asr_tree_card)

    # -------------------------------------------------------------------------
    # Ancestors card

    if ancestor_directory is not None:
        anc_card = create_ancestor_card(anc_dict,
                                        output_directory,
                                        ancestor_directory,
                                        event_color=EVENT_COLOR,
                                        event_name=EVENT_NAME)
        card_stack.append(anc_card)
        

    # Assemble final report html

    container = ["<div class=\"container-lg px-4 gy-4\"><br/>"]
    container.append("".join(card_stack))
    container.append("<br/></div>")
    container_html = "".join(container)

    top, bottom = create_main_html(description=html_description,
                                   title=html_title)

    out = [top,container_html,bottom]

    f = open(os.path.join(output_directory,"index.html"),"w")
    f.write("\n".join(out))
    f.close()
    
    if create_zip_file:
        shutil.make_archive(f"{output_directory}","zip",output_directory)



def pipeline_report(pipeline_directory,
                    output_directory,
                    html_description="Topiary ancestral sequence reconstruction results",
                    html_title="Topiary ASR",
                    overwrite=False,
                    create_zip_file=True):
    """
    Create a sharable html file and directory holding the results of a topiary
    alignment to ancestor (and possible bootstrap reconcile) pipeline.

    Parameters
    ----------
    pipeline_directory : str
        path to directory holding pipeline output
    output_directory : str
        directory in which to write the report
    overwrite : bool, default = False
        whether or not to overwrite an existing report. 
    create_zip_file : bool, default = True
        whether or not to create a zip file of the report directory. This zip
        file will be self-contained and is a way to share the output of the 
        calculation with someone else. 
    """

    # Get the calculation directories
    calc_dirs = _find_directories(pipeline_directory)
    
    model_dir = calc_dirs["model"]
    gene_dirs = calc_dirs["gene"]
    recon_dirs = calc_dirs["reconciled"]

    if model_dir is None and gene_dirs["tree"] is None and recon_dirs["tree"] is None:
        err = f"Pipeline directory '{pipeline_directory}' does not have calculations.\n"
        raise ValueError(err)

    create_output_directory(output_directory,overwrite=overwrite)
    top, bottom = create_main_html(description=html_description,title=html_title)
    
    card_stack = []
    title_html = create_card(html_description,pipeline_directory,"h3")
    card_stack.append(f"{title_html}<br/>")

    tree_stack = []
    for some_dir in [model_dir,gene_dirs["tree"],recon_dirs["tree"]]:
        if some_dir is not None:
            sv = Supervisor(some_dir)
            sv.df.to_csv(os.path.join(output_directory,"dataframe.csv"),index=False)
            input_html = create_input_card(sv)
            tree_stack.append(input_html)
            break

    if gene_dirs["tree"] is not None:

        out_dir = f"{output_directory}/gene-tree/"
        tree_report(tree_directory=gene_dirs["tree"],
                    ancestor_directory=gene_dirs["anc"],
                    output_directory=out_dir,
                    html_description="Ancestors reconstructed on gene tree",
                    html_title="Topiary gene tree ASR",
                    overwrite=overwrite,
                    create_zip_file=False)

        html = f"<br/><h4><a href=\"gene-tree/index.html\">Gene Tree Ancestors</a></h4>"

        help_html = create_info_modal(modal_text=gene_tree_help_text,
                                      modal_title="Gene tree ancestors",
                                      extra_button_class="text-end")
        help_html = 7*"<br/>" + help_html
        html = f"{html}{help_html}"

        tree_stack.append(create_card(html))

    if recon_dirs["tree"] is not None:

        out_dir = f"{output_directory}/reconciled-tree/"
        tree_report(tree_directory=recon_dirs["tree"],
                    ancestor_directory=recon_dirs["anc"],
                    output_directory=out_dir,
                    html_description="Ancestors reconstructed on reconciled gene/species tree",
                    html_title="Topiary reconciled ASR",
                    overwrite=overwrite,
                    create_zip_file=False)

        html = f"<br/><h4 style=\"vertical-align:middle;\" class=\"align-middle\"><a href=\"reconciled-tree/index.html\">Reconciled Tree Ancestors</a></h4>"

        help_html = create_info_modal(modal_text=reconciled_tree_help_text,
                                      modal_title="Reconciled tree ancestors",
                                      extra_button_class="text-end")
        help_html = 7*"<br/>" + help_html
        html = f"{html}{help_html}"

        tree_stack.append(create_card(html))


    if len(tree_stack) > 0:
        card_stack.append(create_row(tree_stack))

    if model_dir is not None:
        model_html = create_model_card(Supervisor(model_dir),output_directory)
        
        container = ["<div class=\"container-lg px-4 gy-4\">"]
        container.append(model_html)
        container.append("</div>")
        container_html = "".join(container)

        #card_stack.append(container_html)
        card_stack.append(model_html)


    container = ["<div class=\"container-lg px-4 gy-4\"><br/>"]
    container.append("".join(card_stack))
    container.append("</br></div>")
    container_html = "".join(container)

    top, bottom = create_main_html(description=html_description,
                                   title=html_title)

    out = [top,container_html,bottom]

    f = open(os.path.join(output_directory,"index.html"),"w")
    f.write("\n".join(out))
    f.close()

    if create_zip_file:
        shutil.make_archive(f"{output_directory}","zip",output_directory)