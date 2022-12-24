
import topiary
from topiary.reports.elements import df_to_table
from topiary.reports.elements import create_card
from topiary.reports.elements import create_icon_row

import numpy as np
import pandas as pd

help_text = \
"""
"""

def create_input_card(supervisor,
                      p_column=None):
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


