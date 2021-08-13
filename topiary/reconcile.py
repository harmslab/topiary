

def annotate_tree_with_calls(df,tree,work_on_copy=True):
    """
    Annotate the leaves of an ete3 tree with information extracted from a
    topiary dataframe.
    """

    # Copy tree -- do not operate on input tree directly
    tree = util.load_tree(tree)

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

def create_generax(df,gene_tree):
    """
    Gene tree, expected to have uid taxon names.
    """

    # Only look at sequences flagged to keep
    df = df.loc[df.keep].copy() # copy is to avoid assign-to-copy warning

    # Annotate gene tree with uid, ott, etc. and arbitrarily resolve any
    # polytomies.
    gene_tree = annotate_tree_with_calls(df,gene_tree)
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
    dp_tree = topiary.get_species_tree(df)
    dp_tree_str = dp_tree.as_string(schema="newick")
    dp_tree_str = re.sub("\'","",dp_tree_str)
    species_tree = ete3.Tree(dp_tree_str)

    # Get species tree that has only the species from the gene tree
    species_tree_to_write = species_tree.copy(method="deepcopy")
    species_tree_to_write.prune(list(seen_in_gene_tree))

    # Resolve polytomies and make sure all branch lenghts/supports have values
    species_tree.resolve_polytomy()
    for n in species_tree_to_write.traverse():
        if n.dist != 1:
            n.dist = 1
        if n.support != 1:
            n.support = 1

    f = open("generax-test/alignment.fasta","w")
    f.write("".join(fasta_out))
    f.close()

    gene_tree_to_write.write(outfile="generax-test/gene_tree.newick",format=5)
    species_tree_to_write.write(outfile="generax-test/species_tree.newick")

    f = open("generax-test/mapping.link","w")
    for k in link_dict:
        f.write(f"{k}:")
        f.write(";".join(link_dict[k]))
        f.write("\n")
    f.close()
