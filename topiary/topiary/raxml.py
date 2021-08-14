__description__ = \
"""
Wrap key features of raxml for asr work.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2021-07-22"

# raxml binary to use it not specified by user
RAXML_BINARY = "raxml-ng.dev"

# -----------------------------------------------------------------------------
# Module import
# -----------------------------------------------------------------------------

import pastml.acr
import ete3
from ete3 import Tree

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec

import subprocess, os, glob, re, sys, time, random, string
import shutil, multiprocessing

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
# raxml interaction functions
# -----------------------------------------------------------------------------

def _parse_raxml_info_for_aic(info_file):
    """
    Open a log file from a tree evaluation run and get log likelihood,
    number of fit parameters, and various AIC scores.
    """

    # Open the file and read lines
    out = {}
    with open(info_file,'r') as f:
        for line in f:

            # Look for likelihood
            if re.search("Final LogLikelihood:",line):
                out["L"] = float(line.strip().split(":")[-1])
                continue

            # Look for number of parameters
            if re.search("Free parameters",line):
                out["N"] = int(line.strip().split(":")[-1])
                continue

            if re.search("AIC score",line):
                cols = line.split("/")
                for c in cols:

                    value = float(c.strip().split(":")[-1])
                    key = c.strip().split(":")[0].split()[0].strip()
                    out[key] = value

    # Return L, N, and AIC
    return out

def _gen_seed():
    """
    Generate a random string of 10 integers and return as a string for passing
    to raxml.
    """

    return "".join([f"{random.choice(range(10)):d}" for _ in range(10)])

def _create_new_dir(dir_name=None):
    """
    Create a new directory.

    dir_name: if specified, name the directory this

    returns name of created directory
    """

    # if dir_name is not specified, build it in a stereotyped fashion
    if dir_name is None:
        rand_name = "".join([random.choice(string.ascii_letters)
                              for _ in range(10)])
        dir_base = os.path.split(RAXML_BINARY)[-1]

        dir_name = f"{dir_base}_{rand_name}"

    # If directory already exists, throw error
    if os.path.exists(dir_name):
        err = f"{dir_name} already exists.\n"
        raise FileExistsError(err)

    # Make directory
    os.mkdir(dir_name)

    return dir_name

def _copy_input_file(input_file,
                     dir_name,
                     file_name=None,
                     make_input_dir=True):
    """
    copy an input file into a directory in a stereotyped way.

    If make_input_dir is specified, copy input_file into dir_name/00_input,
    creating 00_input if necessary.  If make_input_dir is not specified,
    copy in the file as dir_name/{input_file}.

    input_file: file to copy in
    dir_name: copy into dir_name
    file_name: what to call file in new directory. If none, use same name.
    make_input_dir: (bool) make input directory 00_input or not.

    returns name of copied file
    """


    if file_name is None:
        file_name = os.path.split(input_file)[-1]
    file_alone = os.path.split(file_name)[-1]

    # If we are putting this into an input subdirectory
    if make_input_dir:
        input_dir = os.path.join(dir_name,"00_input")
        if not os.path.exists(input_dir):
            os.mkdir(input_dir)
        file_alone = os.path.join("00_input",file_alone)

    # If we are not making an input directory, append input_ to front
    #else:
    #    if not file_alone.startswith("input"):
    #        file_alone = f"input_{file_alone}"

    shutil.copy(input_file,os.path.join(dir_name,file_alone))

    return file_alone

def _run_raxml(algorithm=None,
               alignment_file=None,
               tree_file=None,
               model=None,
               dir_name=None,
               seed=None,
               threads=1,
               raxml_binary=RAXML_BINARY,
               other_args=[]):
    """
    Run raxml. Creates a working directory, copies in the relevant files, runs
    there, and then returns to the previous directory.

    algorithm: algorithm to run (--all, --ancestral, etc.)
    alignment_file: alignment file in .phy format (passed via --msa)
    tree_file: tree file in .newick format (passed via --tree)
    model: model in format recognized by --model
    dir_name: If specified, this will be the name of the working directory.
    seed: true/false, int, or str. If true, pass a randomly generated seed to
          raxml. If int or str, use that as the seed. (passed via --seed)
    threads: number of threads to use (passed via --threads)
    raxml_binary: raxml binary to use
    other_args: list of arguments to pass to raxml
    """

    # Create directory in which to do calculation
    dir_name = _create_new_dir(dir_name=dir_name)

    # Copy alignment and tree files into the directory (if specified)
    if alignment_file is not None:
        alignment_file = _copy_input_file(alignment_file,
                                          dir_name,
                                          file_name="alignment",
                                          make_input_dir=False)
    if tree_file is not None:
        tree_file = _copy_input_file(tree_file,
                                     dir_name,
                                     file_name="tree",
                                     make_input_dir=False)

    # Go into working directory
    cwd = os.getcwd()
    os.chdir(dir_name)

    # Build a command list
    cmd = [raxml_binary]

    if algorithm is not None:
        cmd.append(algorithm)

    if alignment_file is not None:
        cmd.extend(["--msa",alignment_file])

    if tree_file is not None:
        cmd.extend(["--tree",tree_file])

    if model is not None:
        cmd.extend(["--model",model])

    # seed argument is overloaded. Interpret based on type
    if seed is not None:
        if type(seed) is bool:
            cmd.extend(["--seed",_gen_seed()])
        elif type(seed) is int:
            cmd.extend(["--seed",f"{seed:d}"])
        elif type(seed) is str:

            try:
                int(seed)
            except ValueError:
                err = f"seed {seed} could not be interpreted as an int\n"
                raise ValueError(err)

            cmd.extend(["--seed",seed])
        else:
            err = "seed must be True/False, int, or string representation of int\n"
            raise ValueError(err)

    cmd.extend(["--threads",f"{threads:d}"])

    # Put on any custom args
    for a in other_args:
        cmd.append(a)

    # Construct command and dump to std out
    full_cmd = " ".join(cmd)
    print(f"Running '{full_cmd}'")
    sys.stdout.flush()

    # Call subprocess with command
    ret = subprocess.run(cmd,stdout=subprocess.PIPE)
    if ret.returncode != 0:
        err = f"ERROR: raxml returned {ret.returncode}\n\n"
        err += "------------------------------------------------------------\n"
        err += " raxml output \n"
        err += "------------------------------------------------------------\n"
        err += "\n\n"

        err += "".join([line for line in ret.stdout.decode()])

        raise RuntimeError(err)

    # Leave working directory
    os.chdir(cwd)

# -----------------------------------------------------------------------------
# Tree operations
# -----------------------------------------------------------------------------

def _generate_parsimony_tree(alignment_file,
                             dir_name="parsimony-tree",
                             threads=1,
                             raxml_binary=RAXML_BINARY):
    """
    Generate a parsimony tree from an alignment.

    alignment_file: alignment file in .phy format
    dir_name: name to give directory
    threads: number of threads to use
    raxml_binary: raxml binary to use
    """

    _run_raxml(algorithm="--start",
               alignment_file=alignment_file,
               dir_name=dir_name,
               seed=True,
               model="LG",
               threads=threads,
               raxml_binary=raxml_binary,
               other_args=["--tree","pars{1}"])


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


def _make_ancestor_summary_trees(avg_pp_dict,
                                 tree_file_with_labels,
                                 tree_file_with_supports=None):
    """
    Make trees summarizng ASR results.

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

    t_out_pp.write(outfile="ancestors_pp.newick",
                   format=2,format_root_node=True)
    t_out_label.write(outfile="ancestors_label.newick",
                      format=3,format_root_node=True)
    t_out_all.write(outfile="ancestors_all.newick",
                    format=3,format_root_node=True)

    if tree_file_with_supports is not None:
        t_out_all.write(outfile="ancestors_support.newick",
                        format=3,format_root_node=True)


def _parse_raxml_anc_output(anc_prob_file,
                            alignment_file,
                            tree_file_with_labels,
                            tree_file_with_supports=None,
                            dir_name="ancestors",
                            alt_cutoff=0.25,
                            plot_width_ratio=5):
    """
    Parse raxml marginal ancestral state reconstruction output and put out in
    human-readable fashion.

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
    dir_name = _create_new_dir(dir_name)

    anc_prob_file = _copy_input_file(anc_prob_file,dir_name)
    alignment_file = _copy_input_file(alignment_file,dir_name)
    tree_file_with_labels = _copy_input_file(tree_file_with_labels,dir_name)
    if tree_file_with_supports is not None:
        tree_file_with_supports = _copy_input_file(tree_file_with_supports,
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
    df = pd.concat(df_list,ignore_index=True)
    df.to_csv(f"ancestors.csv")

    # Create final tree
    _make_ancestor_summary_trees(avg_pp_dict,
                                 tree_file_with_labels,
                                 tree_file_with_supports)

    os.chdir(cwd)


# -----------------------------------------------------------------------------
# Public functions for doing main raxml calculations
# -----------------------------------------------------------------------------

def find_best_model(alignment_file,
                    tree_file=None,
                    model_matrices=["Blosum62","cpREV","Dayhoff","DCMut","DEN",
                                    "FLU","HIVb","HIVw","JTT","JTT-DCMut","LG",
                                    "mtART","mtMAM","mtREV","mtZOA","PMB",
                                    "rtREV","stmtREV","VT","WAG","LG4M","LG4X",
                                    "PROTGTR"],
                    model_rates=["","G8"],
                    model_freqs=["","FC","FO"],
                    model_invariant=["","IO","IC"],
                    output=None,
                    threads=1,
                    raxml_binary=RAXML_BINARY):
    """
    Find the best phylogentic model to use for tree and ancestor reconstruction
    given an alignment and (possibly) a tree.

    alignment_file: alignment file in .phy format
    tree_file: tree file in newick format. If not specified, parsimony tree
               is generated and used
    model_matrices: list of model matrices to check
    model_rates: ways to treat model rates
    model_freqs: ways to treat model freqs.
    output: directory for output. it none, will generate random name
    threads: number of threads to use
    raxml_binary: raxml binary to use
    """

    # Make sure alignment file exists
    if not os.path.exists(alignment_file):
        err = f"alignment file {alignment_file} does not exist\n"
        raise ValueError(err)

    # Create output directory
    if output is None:
        rand = "".join([random.choice(string.ascii_letters) for _ in range(10)])
        output = f"find_best_model_{rand}"

    dir_name = _create_new_dir(dir_name=output)

    # Copy files into input directory
    alignment_file = _copy_input_file(alignment_file,
                                      dir_name,
                                      make_input_dir=True)
    if tree_file is not None:
        tree_file = _copy_input_file(tree_file,
                                     dir_name,
                                     make_input_dir=True)

    # Move into the output directory
    cwd = os.getcwd()
    os.chdir(dir_name)

    # Generate a parsimony tree if not was specified
    if tree_file is None:
        _generate_parsimony_tree(alignment_file,
                                 dir_name="01_make-parsimony-tree",
                                 threads=threads,
                                 raxml_binary=raxml_binary)
        tree_file = "02_parsimony-tree.newick"
        shutil.copy(os.path.join("01_make-parsimony-tree",
                                 "alignment.raxml.startTree"),
                    tree_file)

    # Dictionary to hold stats for each model
    out = {"model":[]}

    seed = _gen_seed()

    # All possible models, dropping rate, freq, invariant for LG4M and LG4X.
    num_models = (len(model_matrices)-2)*len(model_rates)*len(model_freqs)*len(model_invariant) + 2

    # Go over all combos of the requested matrices, rates, and freqs.
    model_counter = 1
    for matrix in model_matrices:
        for rate in model_rates:
            for freq in model_freqs:
                for invariant in model_invariant:

                    # Check for incompatible matrix/freq/rate combos
                    if matrix in ["LG4M","LG4X"]:
                        if rate != "" or freq != "" or invariant != "":
                            continue

                    # Build model string (for example: LG+G8+FC+IO)
                    model = [matrix,rate,freq,invariant]
                    model = [m for m in model if m != ""]
                    model = "+".join(model)

                    # Print model number we're trying
                    print(f"{model} ({model_counter}/{num_models})")
                    model_counter += 1

                    # Optimize branch lengths etc. on the existing tree
                    _run_raxml(algorithm="--evaluate",
                               alignment_file=alignment_file,
                               tree_file=tree_file,
                               model=model,
                               seed=seed,
                               dir_name="tmp",
                               threads=threads,
                               raxml_binary=raxml_binary)

                    # Grab the info file from this run
                    os.chdir("tmp")

                    # Get results from the info file
                    result = _parse_raxml_info_for_aic("alignment.raxml.log")
                    out["model"].append(model)
                    for r in result:
                        try:
                            out[r].append(result[r])
                        except KeyError:
                            out[r] = [result[r]]

                    # Get out of temporary directory and nuke
                    os.chdir("..")
                    shutil.rmtree("tmp")

    # Create a csv file sorted best to worst aicc
    df = pd.DataFrame(out)

    print(df)

    min_aic = np.min(df.AICc)
    df["p"] = np.exp((min_aic - df.AICc)/2)
    indexer = np.argsort(df.p)[::-1]
    df = df.iloc[indexer,:]
    df.to_csv("model-comparison.csv")

    # Get best model
    best_model = df.model.iloc[0]

    # Write model to a file
    f = open("best-model.txt","w")
    f.write(best_model)
    f.close()

    # Print best model to stdout
    print(f"\n\nBest model: {best_model}\nAICc Prob:{df.p.iloc[0]}\n\n")

    # Leave the output directory
    os.chdir(cwd)

def generate_ml_tree(alignment_file,
                     model,
                     tree_file=None,
                     output=None,
                     threads=1,
                     raxml_binary=RAXML_BINARY,
                     write_bs_msa=True):
    """
    Generate maximum likelihood tree with SH supports from an alignment given
    a substitution model.

    alignment_file: alignment in .phy format
    model: model (e.g. LG+G8).
    tree_file: tree_file in newick format. If not specified, a parsimony tree
               will be generated. used as starting point.
    output: name of output directory.
    threads: number of threads to use
    raxml_binary: what raxml binary to use
    write_bs_msa: whether or not to write out all bootstrap alignments
    """

    # Make sure alignment file exists
    if not os.path.exists(alignment_file):
        err = f"alignment file {alignment_file} does not exist\n"
        raise ValueError(err)

    # Create output directory
    if output is None:
        rand = "".join([random.choice(string.ascii_letters) for _ in range(10)])
        output = f"generate_ml_tree_{rand}"

    dir_name = _create_new_dir(dir_name=output)

    # Copy files into input directory
    alignment_file = _copy_input_file(alignment_file,
                                      dir_name,
                                      make_input_dir=True)
    if tree_file is not None:
        tree_file = _copy_input_file(tree_file,
                                     dir_name,
                                     make_input_dir=True)
    # Move into directory
    cwd = os.getcwd()
    os.chdir(dir_name)


    other_args = ["--bs-trees","autoMRE"]
    if write_bs_msa:
        other_args.append("--bs-write-msa")

    # Run raxml to create tree
    _run_raxml(algorithm="--all",
               alignment_file=alignment_file,
               tree_file=tree_file,
               model=model,
               dir_name="01_make-ml-tree",
               seed=True,
               threads=threads,
               raxml_binary=raxml_binary,
               other_args=other_args)
    tree_file = "02_ml-tree.newick"
    shutil.copy("01_make-ml-tree/alignment.raxml.support",tree_file)

    # Leave working directory
    os.chdir(cwd)


def generate_ancestors(alignment_file,
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

    alignment_file: alignment file (phy)
    model: model (e.g. LG+G8).
    tree_file: tree file to use for reconstruction.
    output: name out output directory.
    threads: number of threads to use
    raxml_binary: what raxml binary to use
    alt_cutoff: cutoff to use for altAll

    creates fasta file and csv file with ancestral sequences, set of ancestor
    plots, and a tree with ancestral names and supports
    """

    # Make sure alignment file exists
    if not os.path.exists(alignment_file):
        err = f"alignment file {alignment_file} does not exist\n"
        raise ValueError(err)

    # Make sure tree file exists
    if not os.path.exists(tree_file):
        err = f"tree file {tree_file} does not exist\n"
        raise ValueError(err)

    # Create output directory
    if output is None:
        rand = "".join([random.choice(string.ascii_letters) for _ in range(10)])
        output = f"generate_ancestors_{rand}"

    dir_name = _create_new_dir(dir_name=output)

    # Copy files into input directory. This will put them in 00_input and keep
    # their original filenames so we have some notion of where they came from.
    alignment_file = _copy_input_file(alignment_file,
                                      dir_name,
                                      make_input_dir=True)
    tree_file = _copy_input_file(tree_file,
                                 dir_name,
                                 make_input_dir=True)

    if tree_file_with_supports is not None:
        tree_file_with_supports = _copy_input_file(tree_file_with_supports,
                                                   dir_name,
                                                   make_input_dir=True)

    # Move into working directory
    cwd = os.getcwd()
    os.chdir(dir_name)

    # Do marginal reconstruction on the tree
    _run_raxml(algorithm="--ancestral",
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
    _parse_raxml_anc_output(anc_prob_file,
                            alignment_file,
                            tree_file_with_labels,
                            tree_file_with_supports,
                            dir_name="02_final-ancestors",
                            alt_cutoff=alt_cutoff)

    # Leave working directory
    os.chdir(cwd)
