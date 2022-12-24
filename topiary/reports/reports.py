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

from topiary.reports.cards import create_ancestor_card
from topiary.reports.cards import create_input_card
from topiary.reports.cards import create_model_card

from topiary.reports.quality import check_duplication

import pandas as pd
import numpy as np

import os
import shutil

EVENT_COLOR = {"D":"#64007F","L":"#BAD316","T":"#407E98","S":"#023E55",None:"#000000"}
EVENT_NAME = {"S":"speciation","D":"duplication","L":"loss","T":"transfer",None:'N/A'}

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

    if os.path.isfile(os.path.join(calculation_directory,"run_parameters.json")):
        return calculation_directory, None

    calc_dirs = {"model"     :None,
                 "gene"      :{"anc":[],
                               "tree":[]},
                 "reconciled":{"anc":[],
                               "tree":[]}}
    for c in os.listdir(calculation_directory):

        this_dir = os.path.join(calculation_directory,c)
        if os.path.isdir(this_dir):
            if os.path.exists(os.path.join(this_dir,"run_parameters.json")):
                sv = Supervisor(this_dir)

                # Run has not finished
                if sv.status != "complete":
                    continue

                # creation time
                try:
                    t = sv.run_parameters["creation_time"]
                except KeyError:
                    t = -1                

                if sv.calc_type == "find_best_model":
                    calc_dirs["model"] = this_dir

                # tree type (will be gene, reconciled, or None)
                tree_type = sv.tree_prefix
                if tree_type is None:
                    continue

                # Is an ancestor dir
                if "ancestors" in sv.calc_type:
                    type_key = "anc"
                else:
                    type_key = "tree"

                calc_dirs[tree_type][type_key].append((t,this_dir))

    for tree_type in ["gene","reconciled"]:
        for type_key in ["anc","tree"]:
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
    
    # Get information about ancestors (pp support, evolutionary event, 
    # bs_support, etc.) Stick into anc_dict
    T = load_trees(supervisor.output_dir,prefix=supervisor.tree_prefix)

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
    
    # Build output directory
    create_output_directory(output_directory=output_directory,
                            overwrite=overwrite)

    supervisor.df.to_csv(os.path.join(output_directory,"dataframe.csv"),index=False)

    # Card stack will hold all of the generated html output
    card_stack = []

    # -------------------------------------------------------------------------
    # Title card

    title_html = create_card(html_description,tree_directory,"h3")
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
    # input card

    input_html = create_input_card(supervisor,p_column)
    
    # -------------------------------------------------------------------------
    # Run parameters

    # Figure out if reconciled
    model = supervisor.model
    if supervisor.tree_prefix == "reconciled":
        reconciled = True
    else:
        reconciled = False

    # Figure out if supports have been calculated
    supports = False
    for k in anc_dict:
        if anc_dict[k]["bs_support"] is not None:
            supports = True
        break

    names = ["Evolutionary model",
             "Reconciled gene and species tree",
             "Supports calculated"]

    values = [model,reconciled,supports]

    if ancestor_directory is not None:
        alt_cutoff = Supervisor(ancestor_directory).run_parameters["alt_cutoff"]
        names.append("Ancestor alt-all cutoff")
        values.append(f"{alt_cutoff:.2f}")

    param_df = pd.DataFrame({"name":names,
                             "value":values})
    param_table = df_to_table(param_df,add_header=False,show_row_numbers=False)
    param_html = create_card("Run parameters",card_contents=param_table)

    # -------------------------------------------------------------------------
    # Icon card

    icon_html = create_icon_row([f"ancestor-data.csv",f"ancestors.fasta",f"summary-tree.pdf"],
                                [f"all ancestors csv",f"all ancestors fasta",f"ancestor tree pdf"])

    icon_html = create_card(card_title="Summary files",
                            card_contents=icon_html)

    # Combine
    card_stack.append(create_row([input_html,param_html,icon_html]))

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
        anc_card = create_ancestor_card(anc_dict,
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
    card_stack.append(title_html)

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

        this_html = create_card(f"<h4><a href=\"gene-tree/index.html\">Gene Tree Ancestors</a></h4>")
        tree_stack.append(this_html)

    if recon_dirs["tree"] is not None:

        out_dir = f"{output_directory}/reconciled-tree/"
        tree_report(tree_directory=recon_dirs["tree"],
                    ancestor_directory=recon_dirs["anc"],
                    output_directory=out_dir,
                    html_description="Ancestors reconstructed on reconciled gene/species tree",
                    html_title="Topiary reconciled ASR",
                    overwrite=overwrite,
                    create_zip_file=False)

        this_html = create_card(f"<h4 style=\"vertical-align:middle;\" class=\"align-middle\"><a href=\"reconciled-tree/index.html\">Reconciled Tree Ancestors</a></h4>")
        tree_stack.append(this_html)


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

    if create_zip_file:
        shutil.make_archive(f"{output_directory}","zip",output_directory)