"""
Generate ancestors and various summary outputs.
"""

import topiary

from ._raxml import run_raxml, RAXML_BINARY
from topiary._private.interface import create_new_dir, copy_input_file
from topiary._private.interface import prep_calc, write_run_information
from topiary._private import check

import pastml.acr
import ete3
from ete3 import Tree

import pandas as pd
import numpy as np

import os, re, glob, shutil

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------

# lists of chemically similar amino acids. the same amino acid can occur
# in different lists.
CHEM_SIMILAR = [["A","C","N","Q","S","T"],
                ["A","V"],
                ["H","K","R"],
                ["D","E"],
                ["F","W","Y"],
                ["I","L","M","V"]]

# -----------------------------------------------------------------------------
# Ancestor processing functions
# -----------------------------------------------------------------------------

def _get_ancestral_gaps(alignment_file,tree_file):
    """
    Get ancestral gaps from an alignment and raxml output tree file. Gaps are
    reconstructed by parsimony using the DOWNPASS algorithm as implemented in
    pastml.

    Parameters
    ----------
    alignment_file : str
        phy file used to generate ancestors in RAxML
    tree_file : str
        output tree file with labeled internal nodes

    Returns
    -------
    gap_anc_dict : dict
        dictionary keying internal node names to lists of True (gap), False (no
        gap), and None (gapping unclear) for each site in that ancestor.
    """

    # Read the alignment file
    counter = 0
    leaf_names = []
    with open(alignment_file) as f:
        for line in f:

            # First line
            if counter == 0:
                col = line.split()
                num_taxa = int(col[0])
                num_sites = int(col[1])
                counter += 1

                char_matrix = np.zeros((num_taxa,num_sites),dtype=np.uint8)
                continue

            # Next, blank line
            if line.strip() == "":
                counter += 1
                continue

            # Alternating leaf id and sequence lines
            if counter % 2 == 0:
                leaf_names.append(line.strip())
                counter += 1
            else:
                index = (counter - 3)//2
                char_matrix[index,:] = np.array([c == "-" for c in line.strip()])
                counter += 1

    # Create a data frame where indexes are leaf names and columns are each gap
    out_dict = {}
    for i in range(char_matrix.shape[1]):
        out_dict[f"g{i}"] = char_matrix[:,i]
    gap_df = pd.DataFrame(out_dict)
    gap_df.index = leaf_names

    # Gaps, named as column names
    gap_names = list(gap_df.columns)

    # Load the tree, keeping the internal node names
    tree = Tree(tree_file,format=1)

    # Reconstruct gaps across tree by parsimony
    acr_result = pastml.acr.acr(tree,gap_df,prediction_method="DOWNPASS")

    # Create dictionary keying ancestor name to gap status across sequence
    gap_anc_dict = {}
    for node in tree.traverse('preorder'):
        if node.name in leaf_names:
            continue

        gap_anc_dict[node.name] = []

        # Go through gaps and grab from node feature
        for g in gap_names:
            state = node.__dict__[g]
            if len(state) == 1:
                if state == {0}:
                    state = False
                else:
                    state = True
            else:
                state = None

            gap_anc_dict[node.name].append(state)

    return gap_anc_dict


def _make_ancestor_summary_trees(df,
                                 avg_pp_dict,
                                 tree_file_with_labels,
                                 tree_file_with_supports=None):
    """
    Make trees summarizing ASR results. Creates two or three newick files:
    + ancestors_label.newick. internal names are labeled with ancestor names
    + ancestors_pp.newick. internal names are avg pp for that ancestor
    + ancestors_support.newick. internal names are branch supports

    Parameters
    ----------
    df : pandas.DataFrame
        topiary data frame
    avg_pp_dict : dict
        map ancestor names to avg ancestor posterior probability
    tree_file_with_labels : str
        output from RAxML that has nodes labeled by their ancestor identity.
        Tree Should also have branch lengths.
    tree_file_with_supports : str
        tree file with branch supports (optional)
    """

    # Create label trees
    t_labeled = Tree(tree_file_with_labels,format=1)

    # Create output trees
    t_out_label = Tree(tree_file_with_labels,format=1)
    t_out_pp = Tree(tree_file_with_labels,format=1)

    # Main iterator (over main labeled tree)
    input_label_iterator = t_labeled.traverse("preorder")

    # These sub iterators will get populated with final tree info
    label_iterator = t_out_label.traverse("preorder")
    pp_iterator = t_out_pp.traverse("preorder")

    if tree_file_with_supports is not None:
        t_out_supports = Tree(tree_file_with_supports,format=0)
        support_iterator = t_out_supports.traverse("preorder")

    # Iterate over main iterator
    for input_label_node in input_label_iterator:

        # Update sub itorators
        pp_node = next(pp_iterator)
        label_node = next(label_iterator)

        if tree_file_with_supports is not None:
            support_node = next(support_iterator)

        # See if this is a leaf and make sure the labeles
        is_label_leaf = input_label_node.is_leaf()
        is_pp_leaf = pp_node.is_leaf()

        # If at least one of the support or label nodes is a leaf...
        if is_label_leaf or is_pp_leaf:

            # If both are not leaves, this is bad news
            if not (is_label_leaf and is_pp_leaf):
                err = "Tree order mismatch (node type)\n"
                raise ValueError(err)
            else:

                # If the nodes are both leaves, make sure they have the
                # same name.
                if input_label_node.name != pp_node.name:
                    err = "Tree order mismatch (node name)\n"
                    err += f"{input_label_node.name} != {pp_node.name}\n"
                    raise ValueError(err)

            # If we get here, this is a leaf node. Continue the iteration.
            continue

        anc = input_label_node.name
        anc_name = re.sub("Node","anc",anc)

        label_node.name = anc_name
        pp_node.support = avg_pp_dict[anc_name]

    t_out_pp.write(format=2,format_root_node=True,outfile="ancestors_pp.newick")
    t_out_label.write(format=3,format_root_node=True,outfile="ancestors_label.newick")

    if tree_file_with_supports is not None:
        t_out_supports.write(format=3,format_root_node=True,outfile="ancestors_support.newick")


def _parse_raxml_anc_output(df,
                            anc_prob_file,
                            alignment_file,
                            tree_file_with_labels,
                            tree_file_with_supports=None,
                            dir_name="ancestors",
                            alt_cutoff=0.25,
                            plot_width_ratio=5):
    """
    Parse raxml marginal ancestral state reconstruction output and put out in
    human-readable fashion. Writes fasta file and csv file with ancestors.
    Writes final newick with three trees: ancestor label, posterior probability,
    and branch support. Creates creates summary plots for all ancestors.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary data frame
    anc_prob_file : str
        ancestor posterior probability file as written out by raxml
    alignment_file : str
        phylip alignment used to create ancestors
    tree_file_with_labels : str
        output newick tree file written by raxml
    tree_file_with_supports : str, optional
        newick tree with supports
    dir_name : str, default="ancestors"
        name for output directory
    alt_cutoff : float, default=.25
        cutoff (inclusive) for identifying plausible alternate states
    plot_width_ratio : float,default=5
        ratio of main and histogram plot widths for ancestors
    """

    # Make directory and copy in files
    dir_name = create_new_dir(dir_name)

    anc_prob_file = copy_input_file(anc_prob_file,dir_name)
    alignment_file = copy_input_file(alignment_file,dir_name)
    tree_file_with_labels = copy_input_file(tree_file_with_labels,dir_name)
    if tree_file_with_supports is not None:
        tree_file_with_supports = copy_input_file(tree_file_with_supports,
                                                  dir_name)

    # Move into output directory
    cwd = os.getcwd()
    os.chdir(dir_name)

    # Get gaps, reconstructed by parsimony
    gap_anc_dict = _get_ancestral_gaps(alignment_file,
                                       tree_file_with_labels)

    # Read pp file into a list of ancestors (anc_list) and posterior
    # probabilities for amino acids at each site
    anc_list = []
    anc_all_pp = []
    last_line = ""
    first_line = True
    with open(anc_prob_file) as f:
        for l in f:
            line = l.strip()

            if first_line:
                tmp_column_names = line.split()
                column_names = []
                for c in tmp_column_names:
                    column_names.append(re.sub("p_","",c))

                aa_list = column_names[3:]

                first_line = False
                continue

            col = line.split()

            if len(anc_list) == 0:
                anc_list.append(col[0])
                anc_all_pp.append([])

            if col[0] != anc_list[-1]:
                anc_list.append(col[0])
                anc_all_pp.append([])

            anc_all_pp[-1].append(np.array([float(c) for c in col[3:]]))


    # Create data structures for output
    out = []
    df_list = []
    avg_pp_dict = {}

    # Go through each ancestor
    for i in range(len(anc_list)):

        for_df = {"anc":[],"site":[],"gap":[],
                  "ml_state":[],"ml_pp":[],
                  "alt_state":[],"alt_pp":[],
                  "site_type":[],"entropy":[]}

        anc = anc_list[i]
        anc_name = re.sub("Node","anc",anc)

        # Get ml and next best seq and posterior probability. Entropy holds
        # the shannon entropy for all posterior probabilities across the
        # row.  If entropy is very high, we have very little confidence --
        # likely a gap.
        anc_ml_seq = []
        anc_ml_pp = []
        anc_alt_seq = []
        anc_alt_pp = []
        entropy = []

        # Go through sites
        for j in range(len(anc_all_pp[i])):

            # Get posterior probability vector for site
            pp = anc_all_pp[i][j].copy()

            # Get entropy, ignoring positions with pp = 0
            tmp_pp = pp[pp != 0.0]
            entropy.append(np.sum(-tmp_pp*np.log(tmp_pp)))

            # Get ml posterior probability and sequence
            anc_ml_pp.append(np.max(pp))
            anc_ml_seq.append(aa_list[np.argmax(pp)])

            # Set max to zero to get next best
            pp[np.argmax(pp)] = 0.0

            # Get second best posterior probability and sequence
            anc_alt_pp.append(np.max(pp))
            anc_alt_seq.append(aa_list[np.argmax(pp)])

            # Add gap!
            if gap_anc_dict[anc][j] == True:
                anc_ml_pp[-1] = 1.0
                anc_ml_seq[-1] = "-"
                anc_alt_pp[-1] = 0.0
                anc_alt_seq[-1] = "-"

        # Convert lists read into numpy arrays
        ml_pp = np.array(anc_ml_pp)
        ml_seq = np.array(anc_ml_seq)
        alt_pp = np.array(anc_alt_pp)
        alt_seq = np.array(anc_alt_seq)

        # Create alt all sequence. This sequence is ML at all states except
        # those with pp >= some alternate cutoff (usually 0.25)
        alt_mask = alt_pp >= alt_cutoff
        alt_all_seq = ml_seq.copy()
        alt_all_seq[alt_mask] = alt_seq[alt_mask]
        alt_all_pp = ml_pp.copy()
        alt_all_pp[alt_mask] = alt_pp[alt_mask]

        # Write out information for data frame
        for_df["anc"] = [anc_name for _ in range(len(anc_all_pp[i]))]
        for_df["site"] = [j for j in range(len(anc_all_pp[i]))]
        for_df["gap"] = gap_anc_dict[anc]

        for_df["ml_state"] = list(ml_seq)
        for_df["ml_pp"] = list(ml_pp)
        for_df["alt_state"] = list(alt_seq)
        for_df["alt_pp"] = list(alt_pp)
        for_df["entropy"] = list(entropy)

        # Classify sites according to semi-intelligent criteria. If gapping
        # unclear by parsimony, call as "possible gap."  If it is a gap, call
        # as "gap." If second best pp is above cutoff, call
        # ambiguous.  Look at ml and alt aa to see if this is a chemically
        # similar or dissimilar ambiguity.  Otherwise, call 'good'
        site_type = []
        num_ambig_gaps = 0
        for j in range(len(anc_all_pp[i])):

            if gap_anc_dict[anc][j] is None:
                site_type.append("possible gap")
                num_ambig_gaps += 1

            elif gap_anc_dict[anc][j]:
                site_type.append("gap")

            elif alt_pp[j] >= alt_cutoff:

                found_match = False
                for ml_set in CHEM_SIMILAR:
                    if ml_seq[j] in ml_set:
                        if alt_seq[j] in ml_set:
                            found_match = True
                            continue
                if found_match:
                    site_type.append("ambig_similar")
                else:
                    site_type.append("ambig_dissimilar")
            else:
                site_type.append("good")

        for_df["site_type"] = site_type

        df_list.append(pd.DataFrame(for_df))

        # Calculate average pp, log sum pp, and number of ambiguous sites
        ml_avgPP = np.mean(ml_pp)
        ml_lnPP = np.sum(np.log(ml_pp))
        alt_avgPP = np.mean(alt_all_pp)
        alt_lnPP = np.sum(np.log(alt_all_pp))
        num_ambig = np.sum(alt_mask)

        # Record avg pp
        avg_pp_dict[anc_name] = ml_avgPP

        # Write ML sequence to fasta
        header = [f">{anc_name}",
                  f"avgPP:{ml_avgPP:.3f}",f"lnPP:{ml_lnPP:.3f}",
                  f"num_ambig:{num_ambig}",f"num_ambig_gaps:{num_ambig_gaps}\n"]
        out.append("|".join(header))
        out.append("".join(ml_seq))
        out.append("\n")

        # Write alt sequence to fasta
        header = [f">{anc_name}_altAll",
                  f"avgPP:{alt_avgPP:.3f}",f"lnPP:{alt_lnPP:.3f}",
                  f"num_ambig:{num_ambig}",f"num_ambig_gaps:{num_ambig_gaps}\n"]
        out.append("|".join(header))
        out.append("".join(alt_all_seq))
        out.append("\n")

        # String for subtitle
        anc_data_string = [f"avgPP: {ml_avgPP:.3f}",
                           f"lnPP: {ml_lnPP:.0f}",
                           f"# ambig: {num_ambig}",
                           f"# ambig gaps: {num_ambig_gaps}"]

        # Plot a summary of the ancestor
        topiary.draw.plot_ancestor_data(df_list[-1],
                                        alt_anc_pp=alt_cutoff,
                                        width_ratio=plot_width_ratio,
                                        anc_name=anc_name,
                                        anc_data_string=", ".join(anc_data_string))


    # Write final fasta file
    f = open(f"ancestors.fasta","w")
    f.write("".join(out))
    f.close()

    # Write ancestor df to csv
    anc_df = pd.concat(df_list,ignore_index=True)
    anc_df.to_csv(f"ancestors.csv")

    # Create final tree
    _make_ancestor_summary_trees(df,
                                 avg_pp_dict,
                                 tree_file_with_labels,
                                 tree_file_with_supports)

    os.chdir(cwd)


def generate_ancestors(previous_dir=None,
                       df=None,
                       model=None,
                       tree_file=None,
                       tree_file_with_supports=None,
                       alt_cutoff=0.25,
                       output=None,
                       overwrite=False,
                       num_threads=-1,
                       raxml_binary=RAXML_BINARY):
    """
    Generate ancestors and various summary outputs. Creates fasta file and csv
    file with ancestral sequences, set of ancestor plots, and a tree with
    ancestral names and supports.

    Parameters
    ----------
    previous_dir : str, optional
        directory containing previous calculation. function will grab the the
        csv, model, and tree from the previous run. If this is not specified,
        :code:`df`, :code:`model`, and :code:`tree_file` arguments must be
        specified.
    df : pandas.DataFrame or str, optional
        topiary data frame or csv written out from topiary df. Will override
        dataframe from `previous_dir` if specified.
    model : str, optional
        model (i.e. "LG+G8"). Will override model from `previous_dir`
        if specified.
    tree_file : str
        tree_file in newick format. Will override tree from `previous_dir` if
        specified.
    tree_file_with_supports : str
        tree file with supports for each node in newick format. Will override
        tree_file_with_supports from `previous_dir` if specified.
    alt_cutoff : float, default=0.25
        cutoff to use for altAll calculation. Should be between 0 and 1.
    output : str, optional
        output directory. If not specified, create an output directory
        with form "generate_ancestors_randomletters"
    overwrite : bool, default=False
        whether or not to overwrite existing output
    num_threads : int, default=-1
        number of threads to use. if -1, use all avaialable
    raxml_binary : str, optional
        what raxml binary to use

    Returns
    -------
    plot : toyplot.canvas or None
        if running in jupyter notebook, return toyplot.canvas; otherwise, return
        None.
    """

    result = prep_calc(previous_dir=previous_dir,
                       df=df,
                       model=model,
                       tree_file=tree_file,
                       other_files=[tree_file_with_supports],
                       output=output,
                       overwrite=overwrite,
                       output_base="generate_ancestors")

    df = result["df"]
    csv_file = result["csv_file"]
    model = result["model"]
    tree_file = result["tree_file"]
    alignment_file = result["alignment_file"]
    tree_file_with_supports = result["other_files"][0]
    starting_dir = result["starting_dir"]
    output = result["output"]

    # Remove labels and supports from internal nodes of tree file to make sure
    # compatible with raxml ancestor inference.
    T = ete3.Tree(tree_file)
    T.write(format=4,outfile=tree_file)

    alt_cutoff = check.check_float(alt_cutoff,
                                   "alt_cutoff",
                                   minimum_allowed=0,
                                   maximum_allowed=1)

    # Do reconstruction on the tree
    cmd = run_raxml(algorithm="--ancestral",
                    alignment_file=alignment_file,
                    tree_file=tree_file,
                    model=model,
                    seed=True,
                    dir_name="working_inference",
                    num_threads=num_threads,
                    raxml_binary=raxml_binary)

    anc_prob_file = os.path.join("working_inference",
                                 "alignment.phy.raxml.ancestralProbs")
    tree_file_with_labels = os.path.join("working_inference",
                                         "alignment.phy.raxml.ancestralTree")

    # Parse output and make something human-readable
    _parse_raxml_anc_output(df,
                            anc_prob_file,
                            alignment_file,
                            tree_file_with_labels,
                            tree_file_with_supports,
                            dir_name="working_analysis",
                            alt_cutoff=alt_cutoff)

    outdir = "output"
    os.mkdir(outdir)

    # tree file
    shutil.copy(tree_file,os.path.join(outdir,"tree.newick"))

    # Write run information
    write_run_information(outdir=outdir,
                          df=df,
                          calc_type="ancestors",
                          model=model,
                          cmd=cmd,
                          outgroup=result["outgroup"])

    # Copy ancestor files into an ancestors directory
    files_to_grab = glob.glob(os.path.join("working_analysis","*.*"))
    anc_out = os.path.join(outdir,"ancestors")
    os.mkdir(anc_out)
    for f in files_to_grab:
        shutil.copy(f,os.path.join(anc_out,os.path.split(f)[-1]))

    print(f"\nWrote results to {os.path.abspath(outdir)}\n")

    # Leave working directory
    os.chdir(starting_dir)

    # Create a plot of the tree
    ret = topiary.draw.ancestor_tree(run_dir=output,
                                     output_file=os.path.join(output,
                                                              "output",
                                                              "summary-tree.pdf"))
    return ret
