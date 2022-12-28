"""
Functions for checking the quality of an ancestral inference.
"""

from topiary.reports.elements import df_to_table
from topiary.reports.elements import create_card
from topiary.reports.elements import create_element
from topiary.reports.elements import create_info_modal

import numpy as np
import pandas as pd

duplications_help_text = \
"""
This warning indicates the evolution of your protein may be more complex than 
can be properly treated by the standard phylogenetic models used by RAxML-NG and
GeneRax. <b>Proceed with extreme caution</b>.

The probabilistic models used in ASR are powerful, but do not capture all
possible evolutionary events. One common problem is incomplete lineage sorting
(ILS), where a gene duplicates but exists as several variants in a population
when speciation occurs. Different duplicates are preserved along the
descendant lineages, meaning this cannot be classified as a simple duplication
or speciation event. ILS is a general problem with all ASR methods and is
specifically noted as being outside the scope of GeneRax (the software topiary
uses for gene/species tree reconciliation). 

Another problem is gene fusion, where different parts of a single gene have
different evolutionary histories. The methods used by topiary all assume a
single genetic history for each protein sequence. If we force such a model to
fit a fused alignment, we will likely end up with a nonsensical evolutionary
tree and meaningless ancestral sequences. 

A standard signal for both ILS and gene fusion is high discordance between the
inferred gene and species trees. This manifests as an unexpectedly high
number of duplication and/or transfer events in the reconciled tree. If, for
example, you are studying a protein family where you expect two paralogs, but
you observe 20 duplication events scattered throughout the tree, there is a good
chance that the evolutionary models used for ASR are not appropriate for your
protein family. 

If your protein has more than one domain, one option would be to try to
reconstruct each domain independently. If the discordance disappears, it's good
evidence for a gene fusion event. If the discordance remains, proceed with
extreme caution. 

One way forward in the face of discordance is to compare the sequences--and
functional characteristics--for any ancestors of interest reconstructed using
either the reconciled gene tree and on the gene tree alone. If the results for
ancestors reconstructed on the two trees differ dramatically, one cannot infer
the ancestral sequence with confidence given standard ASR methods. If the results
for the reconstructions on both trees are similar, it suggests whatever features
you are trying to reconstruct are robust to uncertainty in the tree topology. 
"""

def _check_duplication(supervisor,T,p_column):
    """
    Check the duplications observed in a reconciled tree, flagging problems.

    Parameters
    ----------
    supervisor : topiary.Supervisor instance
        Supervisor with a calculation loaded
    T : ete3.Tree
        tree with elements loaded
    p_column : str
        column in supervisor dataframe holding paralog calls
    
    Returns
    -------
    expected_duplications : int
        number of duplications expected
    duplication_count : int
        number of duplications seen
    df : pandas.DataFrame
        dataframe holding duplication information
    """

    this_df = supervisor.df.loc[supervisor.df.keep,:]
    combined_paralogs = list(np.unique(this_df.loc[:,p_column]))
    combined_paralogs = [p.split("|") for p in combined_paralogs]
    paralogs = []
    for c in combined_paralogs:
        paralogs.extend(c)
    paralogs = np.array(np.unique(paralogs))
    paralogs.sort()

    expected_paralogs = paralogs
    expected_duplications = len(expected_paralogs) - 1

    T_dup = T.copy()
    for leaf in T_dup.get_leaves():
        leaf.add_feature("num_duplications",0)

    excess_duplications = {"ancestor":[],
                           "num_descendants":[]}
        
    # Starts at ancestor and moves up
    duplication_count = 0
    for current_node in T_dup.traverse(strategy="levelorder"):
        
        if not current_node.is_leaf():
            
            if current_node.__dict__["event"] == "D":
                duplication_count += 1
                for leaf in current_node.get_leaves():
                    leaf.num_duplications += 1
                
                if leaf.num_duplications > expected_duplications:
                    excess_duplications["ancestor"].append(current_node.anc_label)
                    excess_duplications["num_descendants"].append(len(current_node.get_leaves()))
    
    df = pd.DataFrame(excess_duplications)
    df = df.sort_values(by="num_descendants",ascending=False)

    return expected_duplications, duplication_count, df

def create_duplications_card(supervisor,T,p_column):
    """
    Create a card warning if there are more duplications than expected.
    
    Parameters
    ----------
    supervisor : topiary.Supervisor instance
        supervisor with calculation loaded
    T : ete3.Tree
        tree with events loaded
    p_column : str
        column in supervisor.df that holds the paralogs calls
    
    Returns
    -------
    html : str
        html warning. If no warning is needed, returns and empty string. 
    """

    expect, obs, duplication_df = _check_duplication(supervisor,T,p_column)

    warning_html = ""
    if obs > expect:

        out = []
        out.append(f"We expected {expect} duplications,")
        out.append(f"but observed {obs} duplications.")
        out.append("These extra duplications could be real, but could")
        out.append("also reflect model violation or incomplete lineage")
        out.append("sorting. The unexpected duplications--and number")
        out.append("of affected descendants--are shown below. Please look for")
        out.append("the purple 'D' events in the tree below. See 'Help' for")
        out.append("more information.")

        txt = " ".join(out)
        out = f"<p>{txt}</p><br/>"
        warning_table = df_to_table(duplication_df,show_row_numbers=False)

        s, e = create_element("div",{"class":"overflow-scroll"})
        warning_txt = f"{out}{s}{warning_table}{e}"

        help_html = create_info_modal(modal_text=duplications_help_text,
                                      modal_title="Excess duplications",
                                      extra_button_class="text-end")
        help_html = f"<br/>{help_html}"

        warning_txt = f"{warning_txt}{help_html}"

        warning_html = create_card("<span style=\"color:red;\">Warning</span>",warning_txt,title_tag="h4")

        warning_html = f"{warning_html}<br/>"

    return warning_html
