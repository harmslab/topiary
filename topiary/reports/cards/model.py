
from topiary.reports.elements import create_element
from topiary.reports.elements import df_to_table
from topiary.reports.elements import create_icon_row
from topiary.reports.elements import create_card

import pandas as pd

import shutil
import os

help_text = \
"""
"""

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