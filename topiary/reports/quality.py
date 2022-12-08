"""
Functions for checking the quality of an ancestral inference.
"""

import numpy as np
import pandas as pd

def check_duplication(supervisor,T,p_column):
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