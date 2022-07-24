"""
Sample sequences from a topiary dataframe in a taxonomically informed way.
"""

import topiary
from topiary.quality import score_alignment
from topiary._private import check

import ete3

import numpy as np

import copy

def _prep_species_tree(df,paralog_column):
    """
    Create a species tree with paralogs loaded into paralogs attribute.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe
    paralog_column : str
        column in dataframe containing paralog calls for each sequence.

    Returns
    -------
    df : pandas.DataFrame
        dataframe with organisms who cannot be placed on the current synthetic
        tree set to keep = False
    annotated_species_tree : ete3.Tree
    """

    if paralog_column not in df.columns:
        err = f"\nparalog_column '{paralog_column}' not in dataframe.\n"
        raise ValueError(err)

    # Make sure df has ott
    if "ott" not in df.columns:
        df = topiary.opentree.get_ott_id(df)

    # Get species tree
    species_tree, dropped = topiary.opentree.get_species_tree(df)

    # Drop sequences for species that cannot be resolved on the tree
    df = df.loc[np.logical_not(df.ott.isin(dropped)),:]

    # If everything is dropped, complain
    if len(df) == 0:
        err = "Could not place any species onto a species tree.\n"
        raise ValueError(err)

    all_paralogs = list(np.unique(df.loc[:,paralog_column]))
    paralogs_seen = dict([(p,[]) for p in all_paralogs])
    for leaf in species_tree.get_leaves():
        leaf.paralogs = copy.deepcopy(paralogs_seen)

        this_df = df.loc[df.loc[:,"uid"].isin(leaf.uid),:]
        for i in range(len(this_df)):
            idx = this_df.index[i]
            leaf.paralogs[this_df.loc[idx,paralog_column]].append(this_df.loc[idx,"uid"])

    return df, species_tree

def _even_paralog_budgeting(T,overall_budget):
    """
    Generate as even as possible of a split between paralogs given the
    overall sequence budget.

    Parameters
    ----------
    T : ete3.Tree
        tree with leaves that have the `paralogs` dictionary holding lists of
        uid with paralog names as keys. Assumes all branch lengths are identical.
    overall_budget : int
        number of sequences to keep adding up all all paralogs.

    Returns
    -------
    paralog_budget : dict
        dictionary with number of sequences to keep for each paralog
    """

    # Get all unique paralogs in the tree
    paralogs = []
    for leaf in T.get_leaves():
        paralogs.extend(list(leaf.paralogs.keys()))
    paralogs = list(set(paralogs))

    paralog_counts = dict([(p,0) for p in paralogs])
    for leaf in T.get_leaves():
        for p in paralogs:
            paralog_counts[p] += len(leaf.paralogs[p])

    paralog_budget = {}
    fx = 1.0/(len(paralogs))
    for p in paralog_counts:
        paralog_budget[p] = int(np.round(fx*overall_budget,0))

    # If rounding dropped a count, add it back in
    total_allocated = sum([paralog_budget[p] for p in paralog_budget])
    if total_allocated < overall_budget:
        paralog_budget[paralogs[0]] += 1

    paralog_budget = _finalize_paralog_budget(paralog_budget,paralog_counts)

    return paralog_budget

def _weighted_paralog_budgeting(T,overall_budget):
    """
    Given some overall number of sequences desired, figure out how many to
    assign to each paralog. Decides how to do this in a tree-weighted
    fashion.

    Parameters
    ----------
    T : ete3.Tree
        tree with leaves that have the `paralogs` dictionary holding lists of
        uid with paralog names as keys. Assumes all branch lengths are identical.
    overall_budget : int
        number of sequences to keep adding up all all paralogs.

    Returns
    -------
    paralog_budget : dict
        dictionary holding total number of sequences to keep for each paralog
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
    paralog_counts = dict([(p,0) for p in paralogs])
    paralog_weights = dict([(p,0) for p in paralogs])
    anc = T.get_common_ancestor(T.get_leaves())
    for leaf in T.get_leaves():
        distance = leaf.get_distance(anc)
        weight = np.power(0.5,distance)

        for p in paralogs:
            num_p = len(leaf.paralogs[p])
            paralog_counts[p] += num_p
            paralog_weights[p] += num_p*weight

    # Given our overall budget and the weights on each paralog seen, give
    # each paralog a total budget (paralog_weight/sum(paralog_weights) * overall).
    paralog_budget = {}
    Q = sum([paralog_weights[p] for p in paralog_weights])
    for p in paralog_weights:
        fx = paralog_weights[p]/Q
        paralog_budget[p] = int(np.round(fx*overall_budget,0))

    paralog_budget = _finalize_paralog_budget(paralog_budget,paralog_counts)

    return paralog_budget

def _finalize_paralog_budget(paralog_budget,paralog_counts):
    """
    The budgeting protocols will potentially give a paralog more budget than
    sequences. The loop below will iteratively partition any extra budget to
    other paralogs.

    Parameters
    ----------
    paralog_budget : dict
        planned paralog budget
    paralog_counts : dict
        actual counts for each paralog

    Returns
    -------
    finished_paralogs : dict
        budget redistributed to match actual number of paralog counts
    """

    finished_paralogs = {}
    while True:

        remaining_budget = 0
        for p in list(paralog_budget.keys()):
            if paralog_budget[p] >= paralog_counts[p]:
                remaining_budget += (paralog_budget[p] - paralog_counts[p])
                finished_paralogs[p] = paralog_counts[p]
                paralog_budget.pop(p)

        # We didn't create a remaining budget -- nothing to partition to other
        # paralogs
        if remaining_budget == 0:
            break

        # We finished off all paralog budgets
        if len(paralog_budget) == 0:
            break

        # Re-divide paralog budget among paralogs that have not finished yet
        local_total = sum([paralog_budget[p] for p in paralog_budget])
        budget_to_allocate = local_total + remaining_budget
        for p in paralog_budget:

            # If no budget allocated at this point to any left, assign even
            # weights.
            if local_total == 0:
                fx = 1/len(paralog_budget)
            else:
                fx = paralog_budget[p]/local_total
            paralog_budget[p] = int(np.round(fx*budget_to_allocate,0))

    # Copy in any budget that was not already popped out to finished
    for p in paralog_budget:
        finished_paralogs[p] = paralog_budget[p]

    return finished_paralogs


def _get_sequence_budgets(T,total_budget):
    """
    Figure out the sequence budget for each node in a tree. Starts at the
    ancestor of the tree and then splits sequences as evenly as possible
    among descendants. Assigns a 'budget' attribute to each node indicating
    how many sequences are budgeted for that node and its descendants.

    Parameters
    ----------
    T : ete3.Tree
        tree with `sequences` list associated with each node.
    total_budget : int
        total number of sequences to budget over the whole tree

    Returns
    -------
    T : ete3.Tree
        input tree updated with .budget and .num_seq attributes
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

                current_node.num_seq = sum([len(l.sequences) for l in current_node.get_leaves()])

                if total_budget > current_node.num_seq:
                    total_budget = current_node.num_seq
                current_node.budget = total_budget

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
    """
    Given the annotated budgets on the tree and the number of sequences at the
    tips, create blocks of sequences that must be merged with one another.

    Parameters
    ----------
    T : ete3.Tree
        tree with budgets and num_seq annotated on all nodes

    Returns
    -------
    merge_blocks : list
        list of blocks to merge with the following form
        [(budget_for_merge,list_of_uid_to_merge,ete3_leaf_this_came_from),
         ...]
    """

    # Work on copy -- we're going to edit tree as we go
    T = T.copy()

    merge_blocks = []
    for counter, leaf in enumerate(T.get_leaves()):

        # Already placed into a merge block
        if leaf.budget == -1:
            continue

        # If we have enough budget to take every sequence on this leaf...
        if leaf.num_seq <= leaf.budget:

            # No sequences -- don't do anything
            if leaf.num_seq == 0:
                continue

            # Record that we want to keep all sequences on this leaf. Put into
            # a merge block that that has the same number of output sequences
            # as we ask to merge.
            merge_uid = list(leaf.sequences)
            merge_blocks.append((len(merge_uid),merge_uid,leaf))
            leaf.budget = -1

        # We're taking a subset of the sequences on this leaf
        else:

            # if the budget is greater than zero, merge the sequences on this
            # leaf into a block of size leaf.buddet
            if leaf.budget > 0:
                merge_uid = list(leaf.sequences)
                merge_blocks.append((leaf.budget,merge_uid,leaf))
                leaf.budget = -1

            # If the budget is zero, we need to go down the tree to find
            # a group of leaves to merge
            else:

                for anc in leaf.get_ancestors():

                    # If we hit an ancestor with a non-zero budget, merge all of
                    # its descendants.
                    if anc.budget > 0:

                        budget = anc.budget
                        merge_uid = []

                        # Go through the leaves that descend from this ancestor
                        for m in anc.get_leaves():
                            if m.budget == -1:
                                err = "We hit a descendant leaf of that\n"
                                err += "has already been placed into a merge. This\n"
                                err += "should not be possible for a monophyletic\n"
                                err += "tree!\n\n"
                                raise RuntimeError(err)

                            # Get uid to merge and set budget of this leaf to
                            # -1 -- now merged
                            merge_uid.extend(m.sequences)
                            m.budget = -1

                        merge_blocks.append((budget,merge_uid,anc))

                        # Break out of get_ancestors loop
                        break


    return merge_blocks

def taxonomic_sample(df,
                     paralog_column="recip_paralog",
                     target_seq_number=500,
                     even_paralog_split=True,
                     key_species=[],
                     within_species_redundancy_cutoff=0.95,
                     sparse_column_cutoff=0.95,
                     align_trim=(0.05,0.95),
                     verbose=False):
    """
    Sample sequences from a topiary dataframe in a taxonomically informed way.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe
    paralog_column : str, default="recip_paralog"
        column holding preliminary paralog calls.
    target_seq_number : int, default=500
        number of sequences to aim for in the final dataset
    even_paralog_split : bool, default=True
        try to keep an even number of paralogs. If False, select the relative
        number of each paralog to keep based on its frequency in the dataframe
    key_species : list, default=[]
        list of species (binomial names) that will not have sequences removed
        (unless there are duplicate sequences for that species with identity >
        within_species_redundancy_cutoff).
    within_species_redundancy_cutoff : float, default=0.95
        merge sequences with sequence identity above cutoff when removing
        redundancy of sequences within species.
    sparse_column_cutoff : float, default=0.95
        when checking alignment quality, a column is sparse if it has gaps in
        more than sparse_column_cutoff sequences.
    align_trim : tuple, default=(0.05,0.95)
        when checking alignment quality, do not score the first and last parts
        of the alignment. Interpreted like a slice, but with percentages.
        (0.0,1.0) would not trim; (0.05,0,98) would trim the first 0.05 off the
        front and the last 0.02 off the back.
    verbose : bool, default=False
        output verbosity
    """

    df = check.check_topiary_dataframe(df)

    try:
        df.loc[:,paralog_column]
    except KeyError:
        err = f"\nparalog_column '{paralog_column}' not found in dataframe.\n\n"
        raise ValueError(err)

    target_seq_number = check.check_int(target_seq_number,
                                        "target_seq_number",
                                        minimum_allowed=1)

    even_paralog_split = check.check_bool(even_paralog_split,
                                          "even_paralog_split")

    key_species = check.check_iter(key_species,
                                   "key_species")

    within_species_redundancy_cutoff = check.check_float(within_species_redundancy_cutoff,
                                                         "within_species_redundancy_cutoff",
                                                         minimum_allowed=0.0,
                                                         maximum_allowed=1.0)

    ## sparse_column_cutoff and align_trim passed directly to score_alignment
    # and validated there

    verbose = check.check_bool(verbose,"verbose")

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

        # Score alignment
        this_df = score_alignment(this_df,
                                  alignment_column="alignment",
                                  align_trim=align_trim,
                                  sparse_column_cutoff=sparse_column_cutoff)

        df.loc[paralog_mask,"sparse_run_length"] = this_df["sparse_run_length"]
        df.loc[paralog_mask,"fx_missing_dense"] = this_df["fx_missing_dense"]

    # Construct species tree annotated with paralogs at tips. If there are
    # species we can't resolve on the tree, drop them from the dataframe
    df, species_tree = _prep_species_tree(df,paralog_column=paralog_column)

    print(f"\nInitial sequence quality control: {starting_keep} --> {np.sum(df.keep)} sequences.\n",flush=True)

    # Get the sequence budget for all paralogs
    if even_paralog_split:
        paralog_budget = _even_paralog_budgeting(species_tree,target_seq_number)
    else:
        paralog_budget = _weighted_paralog_budgeting(species_tree,target_seq_number)

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
