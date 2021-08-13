
import topiary

import numpy as np

import os


def _annotate_tree_with_calls(df,tree,work_on_copy=True):
    """
    Annotate the leaves of an ete3 tree with information extracted from a
    topiary dataframe.
    """
    # Copy tree -- do not operate on input tree directly
    tree = topiary.util.load_tree(tree)
    if work_on_copy:
        tree = tree.copy(method="deepcopy")

    # Create dictionaries mapping uid to species, paralog, ott, and call.
    out_dict = {}
    for i in range(len(df)):

        uid = df.uid.iloc[i]
        species = df.species.iloc[i]
        ott = f"{df.ott.iloc[i]:d}"
        paralog = df.paralog.iloc[i]
        call = f"{ott}|{paralog}"

        out_dict[uid] = {"species":species,
                         "paralog":paralog,
                         "ott":ott,
                         "call":call}

    # Go through tree and assign leaves their calls etc.
    for node in tree.get_leaves():
        try:
            for k in out_dict[node.name]:
                node.add_feature(k,out_dict[node.name][k])
        except KeyError:
            err = f"uid {node.name} in tree but not data frame!\n"
            raise ValueError(err)

    return tree

def create_generax(df,gene_tree,model,out_dir="yo"):
    """
    Gene tree, expected to have uid taxon names.
    """

    # Only look at sequences flagged to keep
    df = df.loc[df.keep].copy() # copy is to avoid assign-to-copy warning

    # Annotate gene tree with uid, ott, etc. and arbitrarily resolve any
    # polytomies.
    gene_tree = _annotate_tree_with_calls(df,gene_tree)
    gene_tree.resolve_polytomy()

    uid_in_gene_tree = []
    link_dict = {}
    for l in gene_tree.get_leaves():

        # For generating mapping file, record which uid are associated with what ott
        try:
            link_dict[l.ott].append(l.name)
        except KeyError:
            link_dict[l.ott] = [l.name]

        # Make sure this uid is actually in the df once and only once
        if np.sum(df.uid == l.name) != 1:
            err = f"tree taxon {l.name} either missing or duplicated in data frame\n"
            raise ValueError(err)

        # Record that we saw this uid
        uid_in_gene_tree.append(l.name)

    # Make df only have uid seen (will automatically trim down to only ott
    # of interest)
    mask = np.array([u in uid_in_gene_tree for u in df.uid],dtype=np.bool)
    df = df.loc[mask]

    # Get species tree.
    species_tree = topiary.get_species_tree(df)

    # Resolve polytomies and make sure all branch lenghts/supports have values
    species_tree.resolve_polytomy()
    for n in species_tree.traverse():
        if n.dist != 1:
            n.dist = 1
        if n.support != 1:
            n.support = 1


    # Construct the control file for generax
    control_out = []
    control_out.append("[FAMILIES]")
    control_out.append("- reconcile")
    control_out.append("starting_gene_tree = gene_tree.newick")
    control_out.append("alignment = alignment.phy")
    control_out.append("mapping = mapping.link")
    control_out.append(f"subst_model = {model}")

    # Write out control file
    f = open(os.path.join(out_dir,"control.txt"),"w")
    f.write("\n".join(control_out))
    f.close()

    # Write out .phy file
    topiary.write_phy(df,os.path.join(out_dir,"alignment.phy"),
                      seq_column="alignment")

    # Write out gene tree
    gene_tree.write(outfile=os.path.join(out_dir,"gene_tree.newick"),
                    format=5)

    # Write out species tree
    species_tree.write(outfile=os.path.join(out_dir,"species_tree.newick"))

    # Write out link file
    f = open(os.path.join(out_dir,"mapping.link"),"w")
    for k in link_dict:
        f.write(f"{k}:")
        f.write(";".join(link_dict[k]))
        f.write("\n")
    f.close()
