__description__ = \
"""
Function for sampling from a topiary dataframe in a taxonomically informed
way.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-17"

import topiary
from topiary.quality import score_alignment
from topiary import _arg_processors

import ete3

import numpy as np

import copy

def _prep_species_tree(df,paralog_column):
    """
    Create a species tree with paralogs loaded into paralogs attribute.

    Parameters
    ----------
        df: topiary dataframe
        paralog_column: column in dataframe containing paralog calls for each
                        sequence.

    Return
    ------
        annotated species tree (ete3)
    """

    if paralog_column not in df.columns:
        err = f"\nparalog_column '{paralog_column}' not in dataframe.\n"
        raise ValueError(err)

    # Make sure df has ott
    if "ott" not in df.columns:
        df = topiary.opentree.get_ott_id(df)

    # Get rid of almost identical sequences within each species
    species_tree = topiary.opentree.get_species_tree(df)

    all_paralogs = list(np.unique(df.loc[:,paralog_column]))
    paralogs_seen = dict([(p,[]) for p in all_paralogs])
    for leaf in species_tree.get_leaves():
        leaf.paralogs = copy.deepcopy(paralogs_seen)

        this_df = df.loc[df.loc[:,"uid"].isin(leaf.uid),:]
        for i in range(len(this_df)):
            idx = this_df.index[i]
            leaf.paralogs[this_df.loc[idx,paralog_column]].append(this_df.loc[idx,"uid"])

    return species_tree

def _divide_budget_among_paralogs(T,overall_budget):
    """
    Given some overall number of sequences desired, figure out how many to
    assign to each paralog. Decides how to do this in a tree-weighted
    fashion.

    Parameters
    ----------
        T: ete3 tree with leaves that have the `paralogs` dictionary holding
           lists of uid with paralog names as keys. Assumes all branch
           lengths are identical.
        overall_budget: integer number of sequences to keep adding up all
                        all paralogs.

    Return
    ------
        paralog_budget dictionary holding total number of sequences to keep for
        each paralog
    """

    # Get all unique paralogs in the tree
    paralogs = []
    for leaf in T.get_leaves():
        paralogs.extend(list(leaf.paralogs.keys()))
    paralogs = list(set(paralogs))

    # Figure out a weighted total for each paralog. The weights have two
    # components. Number of paralogs seen (higher weight) and the distance
    # from the ancestor (lower weight). If an organism has two called copies
    # of some paralog and is four splits from the tree ancestor, it will
    # recieve a weight of 2*(0.5^4) = 0.125. If an organism has one called
    # copy of the same paralog, but directly splits off the ancestor, it
    # will have a weight of 1*(0.5^1) = 0.5. This strategy is to ensure we
    # have an adequte number of seqeunces to capture both early duplications
    # and late expansions.
    paralog_totals = dict([(p,0) for p in paralogs])
    anc = T.get_common_ancestor(T.get_leaves())
    for leaf in T.get_leaves():
        distance = leaf.get_distance(anc)
        weight = np.power(0.5,distance)

        for p in paralogs:
            num_p = len(leaf.paralogs[p])
            paralog_totals[p] += num_p*weight

    # Given our overall budget and the weights on each paralog seen, give
    # each paralog a total budget (paralog_weight/sum(paralog_weights) * overall).
    paralog_budget = {}
    all_total = sum([paralog_totals[p] for p in paralog_totals])
    for p in paralog_totals:
        fx = paralog_totals[p]/all_total
        paralog_budget[p] = int(np.round(fx*overall_budget,0))

    return paralog_budget

def _get_sequence_budgets(T,total_budget):
    """
    Figure out the sequence budget for each node in a tree. Starts at the
    ancestor of the tree and then splits sequences as evenly as possible
    among descendants. Assigns a 'budget' attribute to each node indicating
    how many sequences are budgeted for that node and its descendants.

    Parameters
    ----------
        T: ete3 with `sequences` list associated with each node.
        total_budget: total number of sequences to budget over the whole tree

    Return
    ------
        T updated with .budget and .num_seq attributes
    """

    # Starts at ancestor and moves up
    for current_node in T.traverse(strategy="levelorder"):

        # Get ancestor of clade containing sisters to current group
        sisters = current_node.get_sisters()
        num_sister_clades = len(sisters)
        if num_sister_clades == 1:
            sister_node = sisters[0]
        else:
            # This is the root. Record the total budget available and continue
            if num_sister_clades == 0:
                current_node.budget = total_budget
                current_node.num_seq = sum([len(l.sequences) for l in current_node.get_leaves()])
                continue
            else:
                err = "\nInput tree must be rooted and fully resolved (no polytomies)\n\n"
                raise ValueError(err)

        # If we've already visited this split, don't do anything
        if hasattr(current_node,"budget") and hasattr(sister_node,"budget"):
            continue

        # Get available budget from parent of current split
        parent_node = T.get_common_ancestor(current_node,sister_node)
        avail_budget = parent_node.budget

        # Get the number of sequences down the current and sister splits
        current_seq = sum([len(l.sequences) for l in current_node.get_leaves()])
        current_node.num_seq = current_seq

        sister_seq = sum([len(l.sequences) for l in sister_node.get_leaves()])
        sister_node.num_seq = sister_seq

        # Only one sequence in the budget. Only assign if there is only one sequence
        # left. Otherwise, set to 0.
        if avail_budget == 1:
            if current_seq == 1 and sister_seq == 0:
                current_node.budget = 1
                sister_node.budget = 0
            elif current_seq == 0 and sister_seq == 1:
                current_node.budget = 0
                sister_node.budget = 1
            else:
                current_node.budget = 0
                sister_node.budget = 0

        # Try to partition the budget
        else:

            # We have enough budget to accomodate all sequences.
            if current_seq + sister_seq <= avail_budget:

                current_budget = current_seq
                sister_budget = sister_seq

            # We do not have the budget for all sequences.
            else:

                # Amount each would get in an even split
                even_split = avail_budget//2
                remainder = avail_budget % 2

                # If there are fewer current_seq than what they would get with an even
                # split, keep all of them and assign remaining budget to sister.
                if current_seq < even_split:
                    current_budget = current_seq
                    sister_budget = avail_budget - current_budget

                # If there are fewer sister_seq than what they would get with an even
                # split, keep all of them and assign remaining budget to current.
                elif sister_seq < even_split:
                    sister_budget = sister_seq
                    current_budget = avail_budget - sister_budget

                # Both are above what they would get from an even split. Give them both
                # a budget of even_split, with any remainder going to the node with
                # more descendant sequences.
                else:
                    sister_budget = even_split
                    current_budget = even_split

                    if current_seq >= sister_seq:
                        current_budget += remainder
                    else:
                        sister_budget += remainder

            # Assign the allocated budgets to the node
            current_node.budget = current_budget
            sister_node.budget = sister_budget

    return T

def _identify_merge_blocks(T):

    T = T.copy()
    merge_blocks = []
    for counter, leaf in enumerate(T.get_leaves()):

        # Already placed into a merge block
        if leaf.budget == -1:
            continue

        # If we have enough budget to take every sequence
        if leaf.num_seq <= leaf.budget:

            # No sequences -- don't do anything
            if leaf.num_seq == 0:
                continue

            # Record that we want to merge all sequences on this leaf
            merge_uid = list(leaf.sequences)
            merge_blocks.append((len(merge_uid),merge_uid,leaf))
            leaf.budget = -1

        else:

            # If the leaf has some budget, we're taking a subset of the
            # sequences on this leaf.
            if leaf.budget > 0:
                merge_uid = list(leaf.sequences)
                merge_blocks.append((leaf.budget,merge_uid,leaf))
                leaf.budget = -1

            # If the budget is zero, we need to go down the tree to find
            # a group of leaves to merge
            else:
                for anc in leaf.get_ancestors():

                    # If we have a non-zero budget, merge all descendants
                    if anc.budget > 0:

                        budget = anc.budget
                        merge_uid = []
                        for m in anc.get_leaves():
                            if m.budget == -1:
                                err = "We hit a descendant leaf of that\n"
                                err += "has already been placed into a merge. This\n"
                                err += "should not be possible for a monophyletic\n"
                                err += "tree!\n\n"
                                raise RuntimeError(err)
                            merge_uid.extend(m.sequences)
                            m.budget = -1

                        merge_blocks.append((budget,merge_uid,anc))

                        # Break out of get_ancestors loop
                        break

    # Sanity check on merge
    for l in T.get_leaves():
        if l.num_seq > 0:
            assert l.budget == -1

    return merge_blocks

def taxonomic_sample(df,
                     paralog_column="recip_paralog",
                     target_seq_number=500,
                     key_species=[],
                     within_species_redundancy_cutoff=0.99,
                     sparse_column_cutoff=0.95,
                     sparse_run_length_keep_percentile=0.98,
                     fx_missing_dense_cutoff=0.85,
                     align_trim=(0.05,0.95),
                     verbose=False):
    """
    XX

    Parameters
    ----------
        df: topiary dataframe
        paralog_column: column holding preliminary paralog calls
        target_seq_number: number of sequences to put in final dataset
        key_species: list of species (binomial names) that will not have
                     sequences removed (unless there are duplicate sequences
                     for that species with identity > within_species_redundancy_cutoff).
        within_species_redundancy_cutoff: remove duplicate sequences observed
                                          within a species if their identity is
                                          higher than this cutoff
        sparse_run_length_keep_percentile: within each paralog, keep every
                                           sequence that has a max run length
                                           in mostly gaps columns less than this
                                           percentile over the dataset.
        fx_missing_dense_cutoff: drop any sequence with less than this fraction
                                 of the non-sparse columns as non-gap.


        align_trim: do not score the first and last bits of the alignment.
                    Interpreted like a slice, but with percentages. (0.0,1.0)
                    would not trim; (0.05,0,98) would trim the first 0.05 off
                    the front and the last 0.02 off the back.
        verbose: output verbosity

    """

    df = _arg_processors.process_topiary_dataframe(df)

    try:
        df.loc[:,paralog_column]
    except KeyError:
        err = f"\nparalog_column '{paralog_column}' not found in dataframe.\n\n"
        raise ValueError(err)

    target_seq_number = _arg_processors.process_int(target_seq_number,
                                                 "target_seq_number",
                                                 minimum_allowed=1)

    key_species = _arg_processors.process_iter(key_species,
                                               "key_species")

    within_species_redundancy_cutoff = _arg_processors.process_float(within_species_redundancy_cutoff,
                                                                     "within_species_redundancy_cutoff",
                                                                     minimum_allowed=0.0,
                                                                     maximum_allowed=1.0)

    sparse_run_length_keep_percentile = _arg_processors.process_float(sparse_run_length_keep_percentile,
                                                                      "sparse_run_length_keep_percentile",
                                                                      minimum_allowed=0.0,
                                                                      maximum_allowed=1.0)

    fx_missing_dense_cutoff = _arg_processors.process_float(fx_missing_dense_cutoff,
                                                            "fx_missing_dense_cutoff",
                                                            minimum_allowed=0.0,
                                                            maximum_allowed=1.0)

    fx_missing_dense_cutoff = 1 - fx_missing_dense_cutoff

    ## align_trim passed directly to score_alignment and validated there

    verbose = _arg_processors.process_bool(verbose,"verbose")

    # Drop unkept columns
    df = df.loc[df.keep,:]

    # How many sequences we start with
    starting_keep = np.sum(df.keep)
    print(f"Starting with {starting_keep} sequences.\n",flush=True)

    if within_species_redundancy_cutoff < 1.0:
        print("Removing nearly identical sequences within species.\n",flush=True)
        df = topiary.quality.remove_redundancy(df,
                                               cutoff=within_species_redundancy_cutoff,
                                               only_in_species=True)

        if np.sum(df.keep) == 0:
            err = "redundnacy pass removed all sequences!\n"
            raise ValueError(err)

    # Drop unkept columns
    df = df.loc[df.keep,:]

    # Update dataframe with alignment score information
    df["sparse_run_length"] = 0.0
    df["fx_missing_dense"] = 0.0

    # Iterate over all paralogs.
    paralogs = set(df.loc[:,paralog_column])
    for p in paralogs:

        # Create a dataframe with this paralog only
        paralog_mask = df.loc[:,paralog_column] == p
        this_df = df.loc[paralog_mask,:]

        # Align the sequenes within this paralog
        print(f"Performing initial alignment of paralog {p}.\n",flush=True)
        this_df = topiary.run_muscle(this_df,super5=True)

        print(f"\nRemoving poorly aligned sequences.\n",flush=True)

        # Score alignment
        this_df = score_alignment(this_df,
                                  alignment_column="alignment",
                                  align_trim=align_trim,
                                  sparse_column_cutoff=sparse_column_cutoff)

        df.loc[paralog_mask,"sparse_run_length"] = this_df["sparse_run_length"]
        df.loc[paralog_mask,"fx_missing_dense"] = this_df["fx_missing_dense"]

        # Remove hardest-to-align to align sequences from consideration
        tmp = np.array(this_df.sparse_run_length)
        slicer = int(np.round(len(tmp)*sparse_run_length_keep_percentile,0))
        sparse_run_length_cutoff = tmp[np.argsort(tmp)][slicer]

        keep_mask = np.product((this_df.sparse_run_length < sparse_run_length_cutoff,
                                this_df.fx_missing_dense < fx_missing_dense_cutoff),axis=0)

        df.loc[paralog_mask,"keep"] = np.logical_and(df.loc[paralog_mask,"keep"],
                                                     keep_mask)

        if np.sum(keep_mask) == 0:
            err = f"alignment quality pass removed all sequences for paralog {p}!\n"
            raise ValueError(err)


    print(f"\nInitial sequence quality control: {starting_keep} --> {np.sum(df.keep)} sequences.\n",flush=True)

    # Construct species tree annotated with paralogs at tips
    species_tree = _prep_species_tree(df,paralog_column=paralog_column)

    # Get the sequence budget for all paralogs
    paralog_budget = _divide_budget_among_paralogs(species_tree,target_seq_number)

    print("Approximate number of paralogs to keep")
    for p in paralog_budget:
        print(f"    {p}: {paralog_budget[p]}")
    print("",flush=True)

    # Go through each paralog
    uid_to_keep = []
    for p in paralog_budget:

        # Create a tree that has only paralog p in it
        T = species_tree.copy()
        for leaf in T.get_leaves():
            leaf.sequences = tuple(leaf.paralogs[p])
            leaf.name = f"{leaf.name} ({len(leaf.sequences)})"

        # Figure out how to distribute the total budget across the tree
        _get_sequence_budgets(T,paralog_budget[p])

        # Get list of blocks of sequences to merge
        merge_blocks = _identify_merge_blocks(T)

        # Merge each merge block
        for m in merge_blocks:
            budget = m[0]
            uid = m[1]

            this_mask = df.loc[:,"uid"].isin(uid)
            this_df = df.loc[this_mask,:]

            # If a key_species is in the block, take it.
            in_species_mask = this_df.species.isin(key_species)
            if np.sum(in_species_mask) >= budget:
                this_uid = list(this_df.loc[in_species_mask,"uid"])
            else:

                # Sort by best aligner
                a = np.array(this_df.fx_missing_dense/np.sum(this_df.fx_missing_dense))
                b = np.array(this_df.sparse_run_length/np.sum(this_df.sparse_run_length))

                sort_order = np.argsort(a + b)
                this_uid = np.array(this_df.uid.iloc[sort_order])[:budget]

            if verbose:

                species = this_df.loc[this_df.uid.isin(this_uid),"species"]

                print("---------------------------------------------------------")
                if budget == len(uid):
                    print("KEEPING")
                else:
                    print("MERGING")
                print("---------------------------------------------------------")

                print(m[2])
                print()
                print(f"{p}",",".join(species))
                print(flush=True)

            uid_to_keep.extend(this_uid)


    # Update dataframe keep with merge results
    keep_mask = df.loc[:,"uid"].isin(uid_to_keep)
    df.loc[keep_mask,"keep"] = True
    df.loc[np.logical_not(keep_mask),"keep"] = False
    if "always_keep" in df:
        df.loc[df.always_keep,"keep"] = True

    print("Final number of sequences:")
    print(sum(df.keep))
    print()

    return df
