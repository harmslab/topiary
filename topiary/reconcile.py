
def create_generax_output(df,)

fasta_out = []

gene_tree_to_write = gene_tree.copy(method="deepcopy")
gene_tree_to_write.resolve_polytomy()
seen_in_gene_tree = []

link_dict = {}
for l in gene_tree_to_write.get_leaves():

    # For generating mapping file, record which uid are associated with what ott
    try:
        link_dict[l.ott].append(l.name)
    except KeyError:
        link_dict[l.ott] = [l.name]

    # Get sequence
    aligned_seq = df.loc[df.uid == l.name,"alignment"].iloc[0]

    # Clean up alignment, removing non-aa characters
    aligned_seq = re.sub("[^ACDEFGHIKLMNPQRSTVWYZ-]","-",aligned_seq)
    fasta_out.append(f">{l.name}\n{aligned_seq}\n")

    # Record that we saw this ott
    seen_in_gene_tree.append(l.ott)


# Unique list of ott seen
seen_in_gene_tree = list(set(seen_in_gene_tree))

# Create fully resolved species tree
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
