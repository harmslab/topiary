"""
Generate ancestors and various summary outputs.
"""

import topiary

from ._raxml import run_raxml, RAXML_BINARY

from topiary._private import check
from topiary._private.interface import copy_input_file
from topiary._private.interface import create_new_dir
from topiary._private import Supervisor

from topiary.pastml import get_ancestral_gaps

import ete3

import pandas as pd
import numpy as np

import os
import re
import glob
import shutil

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


def _make_ancestor_summary_trees(df,
                                 avg_pp_dict,
                                 tree_file_with_labels):
    """
    Make trees summarizing ASR results. Creates two newick files:

    + tree_anc-label.newick (internal names are labeled with ancestor names)
    + tree_anc-pp.newick (internal names are avg pp for that ancestor)

    Parameters
    ----------
    df : pandas.DataFrame
        topiary data frame
    avg_pp_dict : dict
        map ancestor names to avg ancestor posterior probability
    tree_file_with_labels : str
        output from RAxML that has nodes labeled by their ancestor identity.
        Tree Should also have branch lengths.
    """

    # Create label trees
    t_labeled = ete3.Tree(tree_file_with_labels,format=1)

    # Create output trees
    t_out_label = ete3.Tree(tree_file_with_labels,format=1)
    t_out_pp = ete3.Tree(tree_file_with_labels,format=1)

    # Main iterator (over main labeled tree)
    input_label_iterator = t_labeled.traverse("preorder")

    # These sub iterators will get populated with final tree info
    label_iterator = t_out_label.traverse("preorder")
    pp_iterator = t_out_pp.traverse("preorder")

    # Iterate over main iterator
    for input_label_node in input_label_iterator:

        # Update sub itorators
        pp_node = next(pp_iterator)
        label_node = next(label_iterator)

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

    t_out_pp.write(format=2,format_root_node=True,outfile="tree_anc-pp.newick")
    t_out_label.write(format=3,format_root_node=True,outfile="tree_anc-label.newick")

def _get_bad_columns(phy_file):
    """
    Get columns that are only {"-","X"}. (RAxML drops these; PastML does not.)

    Parameters
    ----------
    phy_file : str
        .phy file used for the analysis. parses assuming sequence names and
        alignments are on separate lines (what topiary uses and raxml dumps out).

    Returns
    -------
    bad_columns : numpy.ndarray
        integer array indexing bad columns in phy file.
    """

    # Load the phy file into a sequence array
    counter = 0
    name = None
    out = {}
    with open(phy_file) as f:
        for line in f:

            # Skip first two lines
            if counter < 2:
                counter += 1
                continue

            # Alternate between name and sequence lines.
            if name is None:
                name = line.strip()
                continue

            out[name] = line.strip()
            name = None

    # Convert to an array
    seq_array = []
    for k in out:
        seq_array.append(np.array(list(out[k])))
    seq_array = np.array(seq_array)

    # Get columns in array that only have {"-","X"}
    bad_values = set({"-","X"})
    bad_columns = []
    for i in range(seq_array.shape[1]):
        if set(seq_array[:,i]).issubset(bad_values):
            bad_columns.append(i)

    return np.array(bad_columns,dtype=int)


def _parse_raxml_anc_output(df,
                            anc_prob_file,
                            alignment_file,
                            tree_file_with_labels,
                            run_directory="ancestors",
                            alt_cutoff=0.25,
                            plot_width_ratio=5):
    """
    Parse raxml marginal ancestral sequence reconstruction output and put out in
    human-readable fashion. Writes fasta file and csv file with ancestors.
    Writes final newick with two trees: ancestor label, posterior probability.
    Creates creates summary plots for all ancestors.

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
    run_directory : str, default="ancestors"
        name for output directory
    alt_cutoff : float, default=.25
        cutoff (inclusive) for identifying plausible alternate states
    plot_width_ratio : float,default=5
        ratio of main and histogram plot widths for ancestors
    """

    shutil.copy(anc_prob_file,run_directory)
    anc_prob_file = copy_input_file(anc_prob_file,run_directory)
    alignment_file = copy_input_file(alignment_file,run_directory)
    tree_file_with_labels = copy_input_file(tree_file_with_labels,run_directory)

    # Move into output directory
    cwd = os.getcwd()
    os.chdir(run_directory)

    # Get gaps, reconstructed by ACR
    gap_anc_dict = get_ancestral_gaps(alignment_file,tree_file_with_labels)

    # Get indexes of columns that RAXML will have dropped
    bad_columns = _get_bad_columns(alignment_file)

    # Read pp file into a list of ancestors (anc_list) and posterior
    # probabilities for amino acids at each site
    anc_list = []
    anc_all_pp = []
    last_line = ""
    first_line = True
    counter = 0
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

            # First ancestor, initialize
            if len(anc_list) == 0:
                counter = 0
                anc_list.append(col[0])
                anc_all_pp.append([])

            # Subsequent ancestors, where the ancestor name changes from what
            # it was before. Initialize
            if col[0] != anc_list[-1]:
                counter = 0
                anc_list.append(col[0])
                anc_all_pp.append([])

            # Add nan values for the posterior probability for any columns that
            # raxml skipped
            while counter in bad_columns:
                anc_all_pp[-1].append(np.nan*np.ones(20,dtype=float))
                counter += 1

            # Record posterior probability from this row
            anc_all_pp[-1].append(np.array([float(c) for c in col[3:]]))
            counter += 1


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
    anc_df.to_csv(f"ancestor-data.csv")

    # Create final tree
    _make_ancestor_summary_trees(df,
                                 avg_pp_dict,
                                 tree_file_with_labels)

    os.chdir(cwd)


def generate_ancestors(prev_calculation=None,
                       df=None,
                       model=None,
                       gene_tree=None,
                       reconciled_tree=None,
                       alt_cutoff=0.25,
                       calc_dir="ancestors",
                       overwrite=False,
                       num_threads=-1,
                       raxml_binary=RAXML_BINARY):
    """
    Generate ancestors and various summary outputs. Creates fasta file and csv
    file with ancestral sequences, set of ancestor plots, and a tree with
    ancestral names and posterior probabilities. Note: this will always
    reconstruct on the reconciled tree if present (either passed in via
    prev_calculation or via the reconciled_tree argument).

    Parameters
    ----------
    prev_calculation : str or Supervisor, optional
        previously completed calculation. Should either be a directory
        containing the calculation (e.g. the directory with run_parameters.json,
        input, working, output) or a Supervisor instance with a calculation
        loaded. Function will load dataframe, model, gene_tree, and
        reconciled_tree from the previous run. If this is not specified, `df`,
        `model`, and `gene_tree` or `reconciled_tree` arguments must be
        specified.
    df : pandas.DataFrame or str, optional
        topiary data frame or csv written out from topiary df. Will override
        dataframe from `prev_calculation` if specified.
    model : str, optional
        model (i.e. "LG+G8"). Will override model from `prev_calculation`
        if specified.
    gene_tree : str or ete3.Tree or dendropy.Tree
        gene_tree. Reconstruct ancestors on this tree. Will override gene_tree
        from `prev_calculation` if specified. Should be newick with only leaf
        names and branch lengths. If this an ete3 or dendropy tree, it will be
        written out with leaf names and branch lengths; all other data will be
        dropped. NOTE: if reconciled_tree is specified OR is present in the
        prev_calculation, the reconciled tree will take precedence over this
        gene tree.
    reconciled_tree : str or ete3.Tree or dendropy.Tree
        reconciled_tree. Reconstruct ancestors on this tree. Will override
        reconciled_tree from `prev_calculation` if specified. Should be newick
        with only leaf names and branch lengths. If this an ete3 or dendropy
        tree, it will be written out with leaf names and branch lengths; all
        other data will be dropped
    alt_cutoff : float, default=0.25
        cutoff to use for altAll calculation. Should be between 0 and 1.
    calc_dir : str, default="ancestors"
        calculation directory.
    overwrite : bool, default=False
        whether or not to overwrite existing calc_dir
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

    # Load in previous calculation. Three possibilities here: prev_calculation
    # is a supervisor (just use it); prev_calculation is a directory (create a
    # supervisor from it); prev_calculation is None (create an empty supervisor).
    if isinstance(prev_calculation,Supervisor):
        supervisor = prev_calculation
    else:
        supervisor = Supervisor(calc_dir=prev_calculation)

    supervisor.create_calc_dir(calc_dir=calc_dir,
                               calc_type="ancestors",
                               overwrite=overwrite,
                               df=df,
                               gene_tree=gene_tree,
                               reconciled_tree=reconciled_tree,
                               model=model)

    alt_cutoff = check.check_float(alt_cutoff,
                                   "alt_cutoff",
                                   minimum_allowed=0,
                                   maximum_allowed=1)

    supervisor.update("alt_cutoff",alt_cutoff)
    supervisor.check_required(required_values=["model","alt_cutoff"],
                              required_files=["alignment.phy","dataframe.csv"])

    # Figure out which tree to use for the reconstruction
    try:
        supervisor.check_required(required_files=["reconciled-tree.newick"])
        tree_prefix = "reconciled"
        reconstruct_tree = supervisor.reconciled_tree
        supervisor.update("calc_type","reconcile_ancestors")
    except FileNotFoundError:
        try:
            supervisor.check_required(required_files=["gene-tree.newick"])
            tree_prefix = "gene"
            reconstruct_tree = supervisor.gene_tree
            supervisor.update("calc_type","ml_ancestors")
        except FileNotFoundError:
            err = "\nancestral reconstruction requires an existing reconciled\n"
            err += "or gene tree.\n\n"
            raise ValueError(err)

    os.chdir(supervisor.working_dir)

    print("Reconstructing ancestral states.\n",flush=True)

    # Do reconstruction on the tree
    cmd = run_raxml(run_directory="00_inference",
                    algorithm="--ancestral",
                    alignment_file=supervisor.alignment,
                    tree_file=reconstruct_tree,
                    model=supervisor.model,
                    seed=supervisor.seed,
                    log_to_stdout=False,
                    suppress_output=True,
                    supervisor=supervisor,
                    num_threads=num_threads,
                    raxml_binary=raxml_binary)

    anc_prob_file = os.path.join("00_inference",
                                 "alignment.phy.raxml.ancestralProbs")
    tree_file_with_labels = os.path.join("00_inference",
                                         "alignment.phy.raxml.ancestralTree")

    # Parse output and make something human-readable. Do this calculation in
    # calc_dir/working/analysis
    analysis_dir = os.path.join("01_parse")
    os.mkdir(analysis_dir)
    supervisor.event("parsing ancestral output",
                     alt_cutoff=supervisor.run_parameters["alt_cutoff"])

    _parse_raxml_anc_output(df,
                            anc_prob_file,
                            supervisor.alignment,
                            tree_file_with_labels,
                            run_directory=analysis_dir,
                            alt_cutoff=supervisor.run_parameters["alt_cutoff"])


    # Copy ancestors with labels and posterior probabilities
    supervisor.stash(os.path.join("01_parse","tree_anc-label.newick"),
                     target_name=f"{tree_prefix}-tree_anc-label.newick")
    supervisor.stash(os.path.join("01_parse","tree_anc-pp.newick"),
                     target_name=f"{tree_prefix}-tree_anc-pp.newick")

    # Copy ancestor files into an ancestors directory
    files_to_grab = glob.glob(os.path.join("01_parse","*.*"))
    anc_out = os.path.join(supervisor.output_dir,f"{tree_prefix}-tree_ancestors")
    os.mkdir(anc_out)
    for f in files_to_grab:
        shutil.copy(f,os.path.join(anc_out,os.path.basename(f)))

    # Close off
    os.chdir(supervisor.starting_dir)
    return supervisor.finalize(successful=True,plot_if_success=True)
