"""
Reciprocal blast sequences.
"""

import topiary
from topiary._private import check
from .ncbi import ncbi_blast
from .local import local_blast

from topiary.ncbi import parse_ncbi_line

import pandas as pd
import numpy as np

import re
import sys
import copy
import itertools

def _prepare_for_blast(df,
                       paralog_patterns,
                       local_blast_db,
                       ncbi_blast_db,
                       ignorecase,
                       min_call_prob,
                       partition_temp,
                       drop_combo_fx,
                       use_start_end):
    """
    Check sanity of input parameters and compile patterns. Return a validated
    topiary dataframe, list of sequences, dictionary of patterns to search for,
    and sanitized parameters.

    Parameters
    ----------
    df: pandas.DataFrame
        Will pull sequences from df.sequences. If there are 'start' and 'stop'
        columns in the dataframe, only blast sequences between start/top (for
        example: start = 5, stop = 20 would blast sequence[5:20]. To turn this
        off, set use_start_end = False.
    paralog_patterns : dict
        dictionary with paralogs as keys and patterns to look for (as
        lists of strings or compiled regular expressions) as values
    local_blast_db : str or None
        local database against which to blast
    ncbi_blast_db : str or None
        database on ncbi against which to blast
    ignorecase : bool
        whether to ignore letter case when searching blast
    min_call_prob : float, default=0.95
        hits from all paralogs that yield a regular expression match
        are weighted by their relative bit scores. Each paralog is
        assigned a relative probability. This cutoff is the minimum
        probability the best paralog match must have to result in a
        paralog call. Value should be between 0 and 1 (not inclusive),
        where min_call_prob --> 1 increases the stringency.
    partition_temp : float, default=1,
        when calculating posterior probability of the paralog call, use this for
        weighting: 2^(bit_score/partition_temp). partition_temp should be a
        float > 0. A higher value corresponds to a higher stringency. (The bit
        score difference between the best hit and the rest would have to be
        higher to be significant).
    drop_combo_fx : float
        when deciding whether to call a paralog as a combo (i.e. A/B) versus
        singleton (i.e. A), prefer A if pp_A/pp_AB > drop_combo_fx.
    use_start_end : bool
        whether or not to use start/stop columns in
        dataframe (if present) to slice subset of sequence to blast.

    Return
    ------
    df : pandas.DataFrame
        validated dataframe
    sequence_list : list
        list of sequences to use for blast
    paralog_patterns : dict
        valiated paralog_patterns
    min_call_prob : float
    partition_temp : float
    drop_combo_fx : float
    """

    # Make sure dataframe is a topiary dataframe
    df = check.check_topiary_dataframe(df)

    ignorecase = check.check_bool(ignorecase,"ignorecase")

    paralog_patterns = check.check_paralog_patterns(paralog_patterns,
                                                    ignorecase=ignorecase)
    if len(paralog_patterns) == 0:
        err = "\nparalog_patterns must have at least one entry\n"
        raise ValueError(err)

    # Validate blast database arguments
    if ncbi_blast_db is None and local_blast_db is None:
        err = "\nPlease specificy either ncbi_blast_db OR local_blast_db\n\n"
        raise ValueError(err)

    if ncbi_blast_db is not None and local_blast_db is not None:
        err = "\nPlease specificy either ncbi_blast_db OR\n"
        err += "local_blast_db, but not both.\n\n"
        raise ValueError(err)


    min_call_prob = check.check_float(min_call_prob,
                                      "min_call_prob",
                                      minimum_allowed=0,
                                      maximum_allowed=1,
                                      minimum_inclusive=False,
                                      maximum_inclusive=False)

    partition_temp = check.check_float(partition_temp,
                                       "partition_temp",
                                       minimum_allowed=0,
                                       minimum_inclusive=False)

    drop_combo_fx = check.check_float(drop_combo_fx,
                                      "drop_combo_fx",
                                      minimum_allowed=0,
                                      maximum_allowed=1)

    use_start_end = check.check_bool(use_start_end,"use_start_end")


    # Create list of all sequences in dataframe
    sequence_list = []
    for idx in df.index:

        # Do not recip blast sequences where we already have keep == False
        if not df.loc[idx,"keep"]:
            continue

        # Get sequence
        s = df.loc[idx,"sequence"]

        # Try to get start/end (only blasts subset of sequences). If start and
        # end are not defined in this dataframe, take whole sequence
        try:
            a = df.loc[idx,"start"]
            b = df.loc[idx,"end"]
        except KeyError:
            a = 0
            b = None

        # Ignore start/stop if specified.
        if not use_start_end:
            a = 0
            b = None

        # Record this sequence
        try:
            sequence_list.append(s[a:b])
        except TypeError:
            err = "\nstart and end columns must be integers\n\n"
            raise ValueError(err)

        # Warn if the sequence is empty after slicing
        if sequence_list[-1] == "":
            err = "\nsequence has zero length (after slicing with start/end).\n"
            err += f"row: {df.loc[idx,:]}\n\n"
            raise ValueError(err)

    return df, sequence_list, paralog_patterns, min_call_prob, partition_temp, drop_combo_fx


def _run_blast(sequence_list,
               local_blast_db,
               ncbi_blast_db,
               ncbi_taxid,
               hitlist_size,
               e_value_cutoff,
               gapcosts,
               num_threads,
               keep_blast_xml,
               **kwargs):
    """
    Run blast on sequence_sequence list, returning a list of dataframes -- one
    df of hits for each sequence in sequence_list.

    Parameters
    ----------
    sequence_list : list
        list of sequences to use as blast queries.
    local_blast_db : str
        local database against which to blast
    ncbi_blast_db : str
        database on ncbi against which to blast
    ncbi_taxid : int
        limit search to species specified by taxid for an ncbi search.
    e_value_cutoff : float
        minimum allowable e value for a hit
    gapcosts : tuple
        gap costs (must be length 2 tuple of ints)
    num_threads : int
        number of threads to use for blast search. if -1, use all available.
    keep_blast_xml: whether or not to keep blast xml
    kwargs : dict
        extra keyword arguments are passed directly to biopython
        NcbiblastXXXCommandline (for local blast) or qblast (for remote
        blast). These take precedence over anything specified above
        (hitlist_size, for example).

    Return
    ------
    out : list
        list of dataframes with blast hits, one for each sequence in sequence_list
    """

    # NCBI blast
    if ncbi_blast_db:

        try:
            taxid = kwargs.pop("taxid")
            if ncbi_taxid is None:
                ncbi_taxid = taxid
            else:
                err = "\nplease specify the taxid to use for an ncbi blast search\n"
                err += "using the `ncbi_taxid` keyword argument.\n"
                raise ValueError(err)
        except KeyError:
            pass

        # Warn that NCBI blasting can be slow
        w = "\nBlasting against the NCBI database can be slow/unstable. Consider\n"
        w += "creating a local BLAST database for your reciprocal BLAST needs.\n"
        print(w,flush=True)

        hit_dfs = ncbi_blast(sequence_list,
                             db=ncbi_blast_db,
                             taxid=ncbi_taxid,
                             blast_program="blastp",
                             hitlist_size=hitlist_size,
                             e_value_cutoff=e_value_cutoff,
                             gapcosts=gapcosts,
                             num_threads=num_threads,
                             keep_blast_xml=keep_blast_xml,
                             **kwargs)

    # Local blast
    else:
        hit_dfs = local_blast(sequence_list,
                              db=local_blast_db,
                              blast_program="blastp",
                              hitlist_size=hitlist_size,
                              e_value_cutoff=e_value_cutoff,
                              gapcosts=gapcosts,
                              num_threads=num_threads,
                              keep_blast_xml=keep_blast_xml,
                              **kwargs)

    return hit_dfs


def _calc_hit_post_prob(hits,paralog_patterns,partition_temp):
    """
    Calculate the posterior probability for paralog matches.

    Parameters
    ----------
    hits : pandas.DataFrame
        dataframe with blast hits
    paralog_patterns : dict
        dictionary with paralogs as keys and patterns to look for (as
        compiled regular expressions) as values
    partition_temp : float
        partition temperature

    Returns
    -------
    paralogs : numpy.ndarray
        paralogs, sorted alphabetically/numerically
    posterior_prob : numpy.ndarray
        array of posterior probabilities for each paralog
    pattern_masks : list
        list of numpy arrays (bool), one array for each paralog. The True/False
        in each array records whether each row in the dataframe was matched by
        that paralog.
    """

    # Get weighted bits for all hits in the dataframe
    all_weights = np.power(2,hits.loc[:,"bits"]/partition_temp)
    Q = np.sum(all_weights)
    initial_partition_temp = partition_temp
    while np.isinf(Q):
        partition_temp = partition_temp*2
        all_weights = np.power(2,hits.loc[:,"bits"]/partition_temp)
        Q = np.sum(all_weights)

    # if partition_temp != initial_partition_temp:
    #     print(f"Adjusted partition_temp from {initial_partition_temp:.3e} to ")
    #     print(f"{partition_temp:.3e} to avoid a numerical overflow.\n")

    # Go through all patterns
    pattern_masks = []
    posterior_prob = []
    paralogs = list(paralog_patterns.keys())
    paralogs.sort()
    for p in paralogs:

        # Make a mask for whether this paralog_patterns hits or not along each
        # description.
        mask = []
        for i in range(len(hits)):
            description = hits.loc[hits.index[i],"hit_def"]
            if paralog_patterns[p].search(description):
                mask.append(True)
            else:
                mask.append(False)

        # Final mask is value keyed to paralog name
        pattern_masks.append(np.array(mask,dtype=bool))

        # Posterior probability is vector with sums of weights for all
        # hits from this paralog
        posterior_prob.append(np.sum(all_weights[mask]))

    # Weight posterior probability all hits from each paralog to all other
    # hits.
    posterior_prob = np.array(posterior_prob)
    posterior_prob = posterior_prob/(np.sum(all_weights))
    paralogs = np.array(paralogs)

    return paralogs, posterior_prob, pattern_masks



def _make_recip_blast_calls(df,
                            hit_dfs,
                            paralog_patterns,
                            min_call_prob,
                            partition_temp,
                            drop_combo_fx,
                            ncbi_blast_db):
    """
    Make paralog calls given blast output and list of patterns.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe with query sequences
    hit_dfs : list
        list of dataframes returned by blast
    paralog_patterns : dict
        dictionary with paralogs as keys and lists of patterns or compiled
        pattern regular regular expressions as values
    min_call_prob : float
        hits from all paralogs that yield a regular expression match
        are weighted by their relative bit scores. Each paralog is
        assigned a relative probability. This cutoff is the minimum
        probability the best paralog match must have to result in a
        paralog call. Value should be between 0 and 1 (not inclusive),
        where min_call_prob --> 1 increases the stringency.
    partition_temp : float
        when calculating posterior probability of the best versus next-best
        bit score, use this for weighting: 2^(bit_score/partition_temp)
    drop_combo_fx : float
        when deciding whether to call a paralog as a combo (i.e. A/B) versus
        singleton (i.e. A), prefer A if pp_A/pp_AB > drop_combo_fx.
    ncbi_blast_db : str
        database used for ncbi blast. (Used to determine if this should be
        parsed as ncbi vs. local blast inputs)

    Return
    ------
    df : pandas.DataFrame
        dataframe with recip blast calls
    """

    # hit_dfs is a list of dataframes, each of which has the blast hits for
    # each sequence in the input topiary dataframe.

    # Go through hits from each sequence
    results = {"recip_found_paralog":[],
               "recip_hit":[],
               "recip_paralog":[],
               "recip_prob_match":[],
               "recip_bit_score":[]}

    for _, hits in enumerate(hit_dfs):

        # No recip blast hits at all for this sequence
        if len(hits) == 0:
            results["recip_found_paralog"].append(False)
            results["recip_hit"].append(pd.NA)
            results["recip_paralog"].append(pd.NA)
            results["recip_prob_match"].append(np.nan)
            results["recip_bit_score"].append(np.nan)
            continue

        paralogs, posterior_prob, pattern_masks = _calc_hit_post_prob(hits,
                                                                      paralog_patterns,
                                                                      partition_temp)

        # Get all possible combinations of the paralogs
        combinations = []
        indexes = range(len(paralogs))
        for i in range(1,len(indexes)+1):
            combinations.extend(list(itertools.combinations(indexes,i)))

        # Get the total posterior probability each combination
        combo_pp = []
        for i, c in enumerate(combinations):
            idx = np.array(c,dtype=int)
            pp_to_add = posterior_prob[idx]
            combo_pp.append((np.sum(pp_to_add),i))

        # Sort by total posterior probability from best to worst. Combos with
        # more paralogs will always be higher, so this will effectively be
        # ordered from N, N-1, N-2 etc. paralogs included in the combo. This
        # checks to see if a paralog combo with fewer paralogs is basically
        # just as good as a paralog with more.
        combo_pp.sort(reverse=True)
        best_combo_index = 0
        for i in range(1,len(combo_pp)):
            this_pp = combo_pp[best_combo_index][0]
            if not np.isclose(this_pp,0):
                if combo_pp[i][0]/combo_pp[best_combo_index][0] > drop_combo_fx:
                    best_combo_index = i

        # Get best posterior probability, best name of paralog (possibly a
        # combination of more than one paralog separated by "/"), and mask for
        # the overall match (possibly a combination of multiple paralog matches)
        recip_prob_match, best_combo = combo_pp[best_combo_index]
        paralog_indexes = np.array(combinations[best_combo],dtype=int)
        best_paralog = "/".join(paralogs[paralog_indexes])
        best_mask = np.ones(len(pattern_masks[0]),dtype=bool)
        for p in paralog_indexes:
            best_mask = np.logical_or(best_mask,
                                      pattern_masks[p])

        # If the overall posterior probability for the paralog is better than
        # min call prob
        if recip_prob_match > min_call_prob:

            recip_paralog = best_paralog
            recip_found_paralog = True

            # Get best scoring hit for this hit
            recip_hits = hits.loc[best_mask,:]
            best_idx = recip_hits.index[np.argmax(recip_hits.loc[:,"bits"])]
            recip_hit = recip_hits.loc[best_idx,"hit_def"]
            recip_bit_score = recip_hits.loc[best_idx,"bits"]

        else:
            recip_found_paralog = False
            recip_paralog = None

            # Get best scoring hit overall
            best_idx = hits.index[np.argmax(hits.loc[:,"bits"])]
            recip_hit = hits.loc[best_idx,"hit_def"]
            recip_bit_score = hits.loc[best_idx,"bits"]

        # Record results
        results["recip_found_paralog"].append(recip_found_paralog)
        results["recip_hit"].append(recip_hit)
        results["recip_paralog"].append(recip_paralog)
        results["recip_prob_match"].append(recip_prob_match)
        results["recip_bit_score"].append(recip_bit_score)

    # Only load results for values with keep == True
    for k in results:

        # If recip_found_paralog, set to False by default. Overwritten with
        # Trues below
        if k == "recip_found_paralog":
            df[k] = False

        # Otherwise, set to a missing value.
        else:
            df[k] = pd.NA

        # Now set keep with the results from this column
        df.loc[df.keep,k] = results[k]


    # Update keep with recip_found_paralog
    df.loc[:,"keep"] = np.logical_and(df["keep"],df["recip_found_paralog"])

    # If we have sequences set to always_keep, do not drop them even if they
    # do not reciprocal blast properly. Set recip_paralog to their name in case
    # downstream functions use this for naming etc.
    if "always_keep" in df.columns:
        df.loc[:,"keep"] = np.logical_or(df.loc[:,"keep"],df.loc[:,"always_keep"])
        rename = np.logical_and(df.loc[:,"always_keep"],
                                pd.isnull(df.loc[:,"recip_paralog"]))
        df.loc[rename,"recip_paralog"] = df.loc[rename,"name"]

    return df


def recip_blast(df,
                paralog_patterns,
                local_blast_db=None,
                ncbi_blast_db=None,
                ncbi_taxid=None,
                ignorecase=True,
                min_call_prob=0.95,
                partition_temp=1,
                drop_combo_fx=0.9,
                use_start_end=True,
                hitlist_size=10,
                e_value_cutoff=0.01,
                gapcosts=(11,1),
                num_threads=-1,
                keep_blast_xml=False,
                **kwargs):
    """
    Take sequences from a topiary dataframe and do a recip blast analysis
    against an NCBI or local blast database. Looks in blast hits for the
    regular expressions defined in paralog_patterns to call paralog for each
    sequence in df. Returns a copy of the input topiary dataframe with five new
    columns:

    + `recip_found_paralog`: True/False, whether a paralog was found
    + `recip_hit`: string, description for paralog hit (if found) or best hit
      (if no match found)
    + `recip_paralog`: string or None. name of paralog from paralog_patterns
    + `recip_prob_match`: float. probability that this is the correct paralog
      call based on relative evalues of all paralog hits (see Note 2)
    + `recip_bit_score`: float. bit_score for the paralog hit (if found)

    Parameters
    ----------

    df : pandas.DataFrame
        topiary dataframe. Will pull sequences from df.sequences. If there are
        'start' and 'stop' columns in the dataframe, only blast sequences
        between start/top (for example: start = 5, stop = 20 would blast
        sequence[5:20]. To turn this off, set use_start_end = False.
    paralog_patterns : dict
        dictionary with paralogs as keys and lists of patterns or compiled
        pattern regular regular expressions as values. See Note 1 for details
    local_blast_db : str or None, default=None
        Local blast database to use. If None, use an NCBI database. Incompatible
        with ncbi_blast_db.
    ncbi_blast_db : str or None, default=None
        NCBI blast database to use. If None, use a local database. Incompatible
        with local_blast_db.
    ncbi_taxid : str or int or list or None, default=None
        limit ncbi blast search to specified taxid. If None, do not limit
        search. If single str/int, interpret as an integer taxid (i.e. 9606 for
        Homo sapiens). If list of str/int, interpret as multiple taxid (i.e.
        [9606,10090] for Homo sapiens OR Mus musculus).
    ignorecase : bool, default=True
        whether to ignore letter case when searching blast. NOTE: if you pass
        in compiled regular expressions, the flags (such as re.IGNORECASE) on
        the *last* pattern for each paralog will be used for all patterns.
    min_call_prob : float, default=0.95
        hits from all paralogs that yield a regular expression match
        are weighted by their relative bit scores. Each paralog is
        assigned a relative probability. This cutoff is the minimum
        probability the best paralog match must have to result in a
        paralog call. Value should be between 0 and 1 (not inclusive),
        where min_call_prob --> 1 increases the stringency.
    partition_temp : float, default=1,
        when calculating posterior probability of the paralog call, use this for
        weighting: 2^(bit_score/partition_temp). partition_temp should be a
        float > 0. A higher value corresponds to a higher stringency. (The bit
        score difference between the best hit and the rest would have to be
        higher to be significant). This is a minium value: it may be adjusted
        on-the-fly to avoid numerical problems in the calculation. See Note 2.
    drop_combo_fx : float, default=0.90
        when deciding whether to call a paralog as a combination (i.e. A/B)
        versus singleton (i.e. A), prefer A if pp_A/pp_AB > drop_combo_fx.
    use_start_end : bool, default=True
        whether or not to use start/stop columns in dataframe (if present) to
        slice subset of sequence for reciprocal blast.
    hitlist_size : int, default=10
        return only the top hitlist_size hits
    e_value_cutoff : float, default=0.001
        only return hits with e_value better than e_value_cutoff
    gapcosts : tuple, default=(11,1)
        BLAST gapcosts (length 2 tuple of ints)
    num_threads : int, default=-1
        number of threads to use. if -1, use all available. (Multithreading
        rarely speeds up remote BLAST).
    keep_blast_xml : bool, default=False
        whether or not to keep raw blast xml output
    **kwargs : dict, optional
        extra keyword arguments are passed directly to biopython
        blast). These take precedence over anything specified above
        (hitlist_size, for example).

    Returns
    ------
    topiary_dataframe : pandas.DataFrame
        copy of input dataframe with new recip blast columns

    Notes
    -----

    1. :code:`paralog_patterns` are used to match reciprocal blast hits. We
       strongly recommend that you use the patterns generated automatically by
       io.read_seed, as this will find patterns that are unique to each
       paralog. If you want to specify them manually, however, use something
       like:

       .. code-block:: python

            paralog_patterns = {"LY96":["lymphocyte antigen 96","MD-2"],
                                "LY86":["lymphocyte antigen 86","MD-1"]}

       This would mean hits with 'lymphocyte antigen 96' or 'MD-2' will map to
       LY96; hits with 'lymphocyte antigen 86' or 'MD-1' will map to LY86. Hits
       that match neither would not be assigned a paralog. String patterns are
       interpreted literally. The pattern :code:`"[A-Z]"` would be escaped to
       look for :code:`"\[A\-Z\]"`. If you want to use regular expressions in
       your patterns, pass them in as compiled regular expressions. For example,

       .. code-block:: python

            my_pattern = re.compile("[A-Z]")

    2. The posterior probability of a given call is calculated using a weighted
       bit score. For each blast hit, we calculate w = 2^(bit_score/partition_temp).
       We then calculate the posterior probability for each paralog match as
       w_paralog/sum(w_all). See
       MC Frith (2020) Bioinformatics. https://doi.org/10.1093/bioinformatics/btz576
       By increasing the partition_temp parameter, you increase stringency.

    """

    print("Doing reciprocal blast.")

    # Check sanity of input parameters and return a validated topiary dataframe,
    # list of sequences, and compiled set of patterns to search. Note: the
    # values of ncbi_blast_db, local_blast_db, ncbi_taxid, hitlist_size,
    # gapcosts, num_threads, and kwargs are checked by the blast functions
    # themselves.
    out = _prepare_for_blast(df,
                             paralog_patterns,
                             local_blast_db,
                             ncbi_blast_db,
                             ignorecase,
                             min_call_prob,
                             partition_temp,
                             drop_combo_fx,
                             use_start_end)

    df, sequence_list, paralog_patterns, min_call_prob, partition_temp, drop_combo_fx = out

    # Run BLAST on sequence list, returning list of dataframes -- one for each
    # seqeunce in sequence_list
    hit_dfs = _run_blast(sequence_list,
                         local_blast_db,
                         ncbi_blast_db,
                         ncbi_taxid,
                         hitlist_size,
                         e_value_cutoff,
                         gapcosts,
                         num_threads,
                         keep_blast_xml,
                         **kwargs)


    num_keep_at_start = np.sum(df.keep)

    # Make paralog calls given blast output and list of patterns
    out_df = _make_recip_blast_calls(df,
                                     hit_dfs,
                                     paralog_patterns,
                                     min_call_prob,
                                     partition_temp,
                                     drop_combo_fx,
                                     ncbi_blast_db)

    # Write out some summary statistics
    print(f"{np.sum(out_df.keep)} of {num_keep_at_start} sequences passed recip blast.\n")
    print("Found the following numbers of paralogs:")
    for p in paralog_patterns:
        num_of_p = np.sum(out_df.recip_paralog == p)
        print(f"    {p}: {num_of_p}")
    print("",flush=True)

    return out_df
