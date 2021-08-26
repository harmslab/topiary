
def _fix_raxml_tree(raxml_tree,out_file):
    """
    Clean up an raxml [support] newick tree so it is readable by other software.

    raxml_tree: newick file dumped by raxml
    out_file: name of file to write out. (does not check for existance; will
              overwrite)
    """

    # Open raxml tree
    f = open(raxml_tree,"r")
    tree = f.read()
    f.close()

    # Deal with wacky support patterns in raxml output
    support_pattern = re.compile("\):.*?\[.*?\]")
    specific_matches = []
    for x in support_pattern.finditer(tree):
        m = x.group(0)
        support = m.split("[")[1][:-1]
        length = m.split(":")[1].split("[")[0]
        out = f"){support}:{length}"

        p = re.sub("\)","\\\)",m)
        p = re.sub("\[","\\\[",p)
        p = re.sub("\]","\\\]",p)

        specific_matches.append((re.compile(p),out))

    # Actually do substitutions
    for s in specific_matches:
        tree = s[0].sub(s[1],tree)

    # Write output file
    g = open(out_file,"w")
    g.write(tree)
    g.close()
