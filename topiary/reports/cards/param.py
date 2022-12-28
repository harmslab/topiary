from topiary._private import Supervisor

from topiary.reports.elements import create_card
from topiary.reports.elements import df_to_table
from topiary.reports.elements import create_info_modal

import pandas as pd

run_parameters_help_text = \
"""
Shows parameters for this analysis. <b>Evolutionary model</b>: amino acid 
substitution model used. <b>Reconciled gene and species tree</b>: if False, the
tree is a gene tree; if True, the tree is the reconciled gene/species tree. 
<b>Supports calculated</b>: whether or not branch supports have been calculated.
(If running topiary as a pipeline, the "topiary-bootstrap-reconcile" script must
be run to generate branch supports for the reconciled tree). <b>Ancestor 
alt-all cutoff</b>: posterior probability cutoff used for determining the 
sequence of the alt-all ancestor. (This is the dashed line in the ancestor
posterior probability plots).
"""

def create_param_card(supervisor,anc_dict,ancestor_directory):

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
    help_html = create_info_modal(modal_text=run_parameters_help_text,
                                  modal_title="Run parameters",
                                  extra_button_class="text-end")
    help_html = 2*"<br/>" + help_html
    param_table = f"{param_table}{help_html}"

    param_html = create_card("Run parameters",card_contents=param_table,title_tag="h4")

    return param_html