__description__ = \
"""
Generate maximum likelihood ancestors using raxml.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2021-07-22"

import topiary

from ._raxml import create_new_dir, copy_input_file, prep_calc
from ._raxml import run_raxml, RAXML_BINARY

import pastml.acr
import ete3
from ete3 import Tree

import pandas as pd
import numpy as np

import os, re

from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec

# -----------------------------------------------------------------------------
# Configure plotting
# -----------------------------------------------------------------------------

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

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

    alignment_file: phy file used to generate ancestors in RAxML
    tree_file: output tree file with labeled internal nodes

    output: dictionary keying internal node names to lists of True (gap),
            False (no gap), and None (gapping unclear) for each site in that
            ancestor.
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

def _plot_ancestor_data(df_anc,
                        alt_anc_pp,
                        width_ratio,
                        anc_name,
                        anc_data_string):
    """
    Create a summary plot for an ancestor.

    df_anc: ancestral data frame
    alt_anc_pp: cutoff (inclusive) for identifying plausible alternate states
    width_ratio: width ratio for plot
    anc_name: name of ancestor (name of plot, title on graph)
    anc_data_string: data to dump in subtitle
    """

    def _draw_histogram(values,ax,bin_size=0.05,color="gray"):
        """
        Draw a histogram sideways next to main plot.

        values: array to bin
        ax: axis on which to make plot
        bin_size: width of bins (from 0 to 1.0)
        color: color to make histogram bars

        returns maximum number of counts (for constructing xlim later)
        """

        # Create histogram
        counts, bins = np.histogram(values,
                                    bins=np.arange(0,1 + bin_size,bin_size))

        # Draw bars for histogram
        for i in range(len(counts)):
            rect = patches.Rectangle((0,bins[i]),
                                     width=counts[i],
                                     height=(bins[i+1]-bins[i]),
                                     linewidth=1,
                                     edgecolor='black',
                                     facecolor=color,
                                     alpha=0.5)
            ax.add_patch(rect)

        return np.max(counts)

    # Data frames containing unambiguous gaps and other sites
    df_gap = df_anc.loc[df_anc.site_type == "gap",:]
    df_nogap = df_anc.loc[df_anc.site_type != "gap",:]

    # create a figure
    fig = plt.figure()
    fig.set_figheight(4)
    fig.set_figwidth(10)

    # create grid for main and histogram subplots
    spec = gridspec.GridSpec(ncols=2, nrows=1,
                             width_ratios=[width_ratio, 1],
                             wspace=0.01)
    # Greate actual subplots
    ax = [fig.add_subplot(spec[0])]
    ax.append(fig.add_subplot(spec[1],sharey=ax[0]))

    # Plot gaps on main figure. Draw boxes for contiguous gap regions, Code
    # below is too clever by half, but returns contiguous blocks of gaps
    sites = np.array(df_gap.site,dtype=np.uint)
    contiguous = np.split(sites,np.where(np.diff(sites) != 1)[0]+1)

    # Iterate over contiguous blocks
    for c in contiguous:

        # Find gap width.  Minimum gap width is 1.
        width = c[-1] - c[0]
        if width == 0:
            width = 1

        # Draw rectangle for gap block
        rect = patches.Rectangle((c[0],0),width,1,
                                 linewidth=1,
                                 edgecolor='lightgray',
                                 facecolor='lightgray')
        ax[0].add_patch(rect)

    # Create list of ambiguous gaps and draw veritical purple lines at these
    # positions.
    ambig_df = df_anc.loc[df_anc.site_type == "possible gap",:]
    for i in range(len(ambig_df)):
        row = ambig_df.iloc[i]
        ax[0].plot((row.site,row.site),(0,1.05),"--",lw=1,color="purple")

    # Plot ML and alt pp points
    ax[0].plot(df_nogap.site,df_nogap.ml_pp,".",color="black",markersize=8)
    ax[0].plot(df_nogap.site,df_nogap.alt_pp,".",color="red",markersize=8)

    # Plot ML and alt pp lines. Only draw lines over contiguous stretches.
    sites = np.array(df_nogap.site,dtype=np.uint)
    contiguous = np.split(sites,np.where(np.diff(sites) != 1)[0]+1)
    for c in contiguous:
        ax[0].plot(df_anc.site.iloc[c],df_anc.ml_pp.iloc[c],color="black",lw=2)
        ax[0].plot(df_anc.site.iloc[c],df_anc.alt_pp.iloc[c],color="red",lw=2)

    # Plot alt-all cutoff line
    ax[0].plot((np.min(df_anc.site),np.max(df_anc.site)),
               (alt_anc_pp,alt_anc_pp),"--",color="gray")

    # Draw histograms for ml and alt pp in right plot
    max_ml_counts = _draw_histogram(df_nogap.ml_pp,ax[1],color="gray")
    max_alt_counts = _draw_histogram(df_nogap.alt_pp,ax[1],color="red")

    # Figure out biggest value seen in histogram plot
    hist_x_max = np.max((max_ml_counts,max_alt_counts))

    # Plot alt-all cutoff line on histogram
    ax[1].plot((0,hist_x_max),
               (alt_anc_pp,alt_anc_pp),"--",color="gray")

    # Clean up axes for main plot
    ax[0].set_xlabel("alignment site")
    ax[0].set_ylabel("posterior probability")
    ax[0].set_xlim(np.min(df_anc.site),np.max(df_anc.site))
    ax[0].spines["right"].set_visible(False)
    ax[0].spines["top"].set_visible(False)
    ax[0].set_ylim(-0.05,1.1)

    # Clean up axes for histogram plot
    ax[1].set_xlim(0,hist_x_max*1.1)
    ax[1].spines["right"].set_visible(False)
    ax[1].spines["top"].set_visible(False)
    ax[1].spines["bottom"].set_visible(True)
    ax[1].spines["left"].set_visible(False)
    ax[1].axis('off')

    # Main plot title
    fig.suptitle(f"{anc_name}")

    # This bit of wackiness finds where on main plot the midpoint for the
    # whole graph is.
    main_plot_x_length = np.max(df_anc.site) - np.min(df_anc.site)
    total_plot_length = ((1 + width_ratio)/(width_ratio))*main_plot_x_length
    mid_point = total_plot_length/2.0

    # Plot the subtitle on the main plot
    ax[0].text(mid_point,1.08,anc_data_string,ha='center')

    # Save out figure
    fig.savefig(f"{anc_name}.pdf",bbox_inches="tight")


def _make_ancestor_summary_trees(df,
                                 avg_pp_dict,
                                 tree_file_with_labels,
                                 tree_file_with_supports=None):
    """
    Make trees summarizng ASR results.

    df: topiary data frame
    avg_pp_dict: dictionary mapping ancestor names to avg ancestor posterior
                 probability
    tree_file_with_labels: output from RAxML that has nodes labeled by their
                           ancestor identity. Should also have branche lengths.
    tree_file_with_supports: tree file with supports (optional)

    Creates three or four newick files:
        ancestors_label.newick: tree where internal names are labeled with\
                                ancestor names
        ancestors_pp.newick: tree where supports are avg pp for that ancestor
        ancestors_support.newick: tree with supports (optional)
        ancestors_all.newick: tree where internal names are name|pp OR name|pp|support
    """

    # Create label trees
    t_labeled = Tree(tree_file_with_labels,format=1)

    # Create output trees
    t_out_label = Tree(tree_file_with_labels,format=1)
    t_out_pp = Tree(tree_file_with_labels,format=1)
    t_out_all = Tree(tree_file_with_labels,format=1)

    # Main iterator (over main labeled tree)
    input_label_iterator = t_labeled.traverse("preorder")

    # These sub iterators will get populated with final tree info
    label_iterator = t_out_label.traverse("preorder")
    pp_iterator = t_out_pp.traverse("preorder")
    all_iterator = t_out_all.traverse("preorder")

    if tree_file_with_supports is not None:
        t_out_supports = Tree(tree_file_with_supports,format=0)
        support_iterator = t_out_supports.traverse("preorder")

    # Iterate over main iterator
    for input_label_node in input_label_iterator:

        # Update sub itorators
        pp_node = next(pp_iterator)
        label_node = next(label_iterator)
        all_node = next(all_iterator)

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
        try:
            anc_name = f"anc{int(anc):04d}"
        except ValueError:
            anc_name = f"anc{anc}"

        label_node.name = anc_name
        pp_node.support = avg_pp_dict[anc_name]

        combo = [f"{anc_name}",
                 f"{pp_node.support:.2f}"]

        # If support specified
        if tree_file_with_supports is not None:
            combo.append(f"{support_node.support:.0f}")

        all_node.name = "|".join(combo)

    topiary.util.uid_to_pretty(df,
                               t_out_pp.write(format=2,format_root_node=True),
                               out_file="ancestors_pp.newick")

    topiary.util.uid_to_pretty(df,
                               t_out_label.write(format=3,format_root_node=True),
                               out_file="ancestors_label.newick")
    topiary.util.uid_to_pretty(df,
                               t_out_all.write(format=3,format_root_node=True),
                               out_file="ancestors_all.newick")

    if tree_file_with_supports is not None:
        topiary.util.uid_to_pretty(df,
                                   t_out_all.write(format=3,format_root_node=True),
                                   out_file="ancestors_support.newick")


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
    human-readable fashion.

    df: topiary data frame
    anc_prob_file: ancestor posterior probability file as written out by raxml
    alignment_file: phylip alignment used to create ancestors
    tree_file_with_labels: output newick tree file written by raxml
    tree_file_with_supports: newick tree with supports (optional)
    name: name for output directory
    alt_cutoff: cutoff (inclusive) for identifying plausible alternate states
    plot_width_ratio: ratio of main and histogram plot widths for ancestors

    Writes fasta file and csv file with ancestors. Writes final newick with
    three trees: ancestor label, posterior probability, and SH support.  Creates
    creates summary plots for all ancestors.

    returns None
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
        gap_list = gap_anc_dict[anc]

        try:
            anc_name = f"anc{int(anc):04d}"
        except ValueError:
            anc_name = f"anc{anc}"

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
                anc_alt_pp[-1] = 1.0
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
        _plot_ancestor_data(df_list[-1],
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


def generate_ancestors(df,
                       model,
                       tree_file,
                       tree_file_with_supports=None,
                       output=None,
                       threads=1,
                       raxml_binary=RAXML_BINARY,
                       alt_cutoff=0.25,
                       calculate_supports=False):
    """
    Generate ancestors and various summary outputs.

    df: topiary data frame or csv written out from topiary df
    model: model (e.g. LG+G8).
    tree_file: tree file to use for reconstruction.
    output: name out output directory.
    threads: number of threads to use
    raxml_binary: what raxml binary to use
    alt_cutoff: cutoff to use for altAll

    creates fasta file and csv file with ancestral sequences, set of ancestor
    plots, and a tree with ancestral names and supports
    """

    result = prep_calc(df=df,
                       output=output,
                       other_files=[tree_file,tree_file_with_supports],
                       output_base="generate_ancestors")

    df = result["df"]
    csv_file = result["csv_file"]
    alignment_file = result["alignment_file"]
    tree_file = result["other_files"][0]
    tree_file_with_supports = result["other_files"][1]
    starting_dir = result["starting_dir"]

    # Do marginal reconstruction on the tree
    run_raxml(algorithm="--ancestral",
              alignment_file=alignment_file,
              tree_file=tree_file,
              model=model,
              seed=True,
              dir_name="01_calc-marginal-anc",
              threads=threads,
              raxml_binary=raxml_binary)

    anc_prob_file = "01_calc-marginal-anc/alignment.raxml.ancestralProbs"
    tree_file_with_labels = "01_calc-marginal-anc/alignment.raxml.ancestralTree"

    # Parse output and make something human-readable
    _parse_raxml_anc_output(df,
                            anc_prob_file,
                            alignment_file,
                            tree_file_with_labels,
                            tree_file_with_supports,
                            dir_name="02_final-ancestors",
                            alt_cutoff=alt_cutoff)

    # Leave working directory
    os.chdir(starting_dir)
