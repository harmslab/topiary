"""
Reciprocal blast sequences.
"""

import topiary
from topiary._private import check
from .ncbi import ncbi_blast
from .local import local_blast

from topiary.external.ncbi import parse_ncbi_line

import pandas as pd
import numpy as np

import re, sys, copy

def _prepare_for_blast(df,
                       paralog_patterns,
                       local_blast_db,
                       ncbi_blast_db,
                       ignorecase,
                       max_del_best,
                       min_call_prob,
                       use_start_end):
    """
    Check sanity of input parameters and compile patterns. Return a validated
    topiary dataframe, list of sequences, compiled set of patterns to search,
    and sanitized parameters.

    Parameters
    ----------
    df: pandas.DataFrame
        Will pull sequences from df.sequences. If there are 'start' and 'stop'
        columns in the dataframe, only blast sequences between start/top (for
        example: start = 5, stop = 20 would blast sequence[5:20]. To turn this
        off, set use_start_end = False.
    paralog_patterns : dict
        dictionary with paralogs as values and lists of patterns to look for as
        values.
    local_blast_db : str or None
        local database against which to blast
    ncbi_blast_db : str or None
        database on ncbi against which to blast
    ignorecase : bool
        whether to ignore letter case when searching blast
    max_del_best : float
        maximum value for log(e_hit_match) - log(e_hit_best) that
        allows for paralog call. This means the matched recip blast
        hit does not have to be the best hit, but must be within this
        e-value difference of the best hit. A higher number means a
        less stringent cutoff. A value of 0 would require the matched
        hit be the best hit.
    min_call_prob: float
        hits from all paralogs that yield a regular expression match
        are weighted by their relative e-values. Each paralog is
        assigned a relative probability. This cutoff is the minimum
        probability the best paralog match must have to result in a
        paralog call. Value should be between 0 and 1 (not inclusive),
        where min_call_prob --> 1 increases the stringency.
    use_start_end : bool
        whether or not to use start/stop columns in
        dataframe (if present) to slice subset of sequence to blast.

    Return
    ------
    df : pandas.DataFrame
        validated dataframe
    sequence_list : list
        list of sequences to use for blast
    patterns : dict
        valiated paralog patterns
    max_del_best : float
    min_call_prob : float
    """

    # Check type of df
    if type(df) is not pd.DataFrame:
        err = "\ndf should be a topiary dataframe\n\n"
        raise ValueError(err)

    # Check ignorecase
    try:
        ignorecase = bool(int(ignorecase))
    except (ValueError,TypeError):
        err = "\nignorecase must be True or False\n\n"
        raise ValueError(err)

    patterns = check.check_paralog_patterns(paralog_patterns,
                                            ignorecase=ignorecase)
    if len(patterns) == 0:
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

    # Check max_del_best
    try:
        max_del_best = float(max_del_best)
        if max_del_best < 0 or np.isnan(max_del_best):
            raise ValueError

    except (ValueError,TypeError):
        err = "\nmax_del_best must be a number between 0 and positive infinity\n\n"
        raise ValueError(err)

    # Check min_call_prob
    try:
        min_call_prob = float(min_call_prob)
        if min_call_prob <= 0 or min_call_prob >= 1:
            raise ValueError
        if np.isnan(min_call_prob):
            raise ValueError
    except (ValueError,TypeError):
        err = "\nmin_call_prob must be a number between 0 and 1 (not inclusive)\n\n"
        raise ValueError(err)

    # Check use_start_end
    try:
        use_start_end = bool(int(use_start_end))
    except (ValueError,TypeError):
        err = "\nuse_start_end must be True or False\n\n"
        raise ValueError(err)

    # Make sure dataframe is a topiary dataframe
    df = check.check_topiary_dataframe(df)

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

    return df, sequence_list, patterns, max_del_best, min_call_prob


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


def _make_recip_blast_calls(df,
                            hit_dfs,
                            patterns,
                            max_del_best,
                            min_call_prob,
                            ncbi_blast_db):
    """
    Make paralog calls given blast output and list of patterns.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary dataframe with query sequences
    hit_dfs : list
        list of dataframes returned by blast
    patterns : dict
        dictionary with paralogs as values and lists of patterns to look for as
        values.
    max_del_best : float
        maximum value for log(e_hit_match) - log(e_hit_best) that
        allows for paralog call. This means the matched recip blast
        hit does not have to be the best hit, but must be within this
        e-value difference of the best hit. A higher number means a
        less stringent cutoff. A value of 0 would require the matched
        hit be the best hit.
    min_call_prob : float
        hits from all paralogs that yield a regular expression match
        are weighted by their relative e-values. Each paralog is
        assigned a relative probability. This cutoff is the minimum
        probability the best paralog match must have to result in a
        paralog call. Value should be between 0 and 1 (not inclusive),
        where min_call_prob --> 1 increases the stringency.
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
               "recip_del_best":[]}

    for _, hits in enumerate(hit_dfs):

        # No recip blast hits at all for this sequence
        if len(hits) == 0:
            results["recip_found_paralog"].append(False)
            results["recip_hit"].append(pd.NA)
            results["recip_paralog"].append(pd.NA)
            results["recip_prob_match"].append(np.nan)
            results["recip_del_best"].append(np.nan)
            continue

        # Get e value of top hit
        top_e_value = hits.loc[hits.index[0],"e_value"]
        top_hit_def = hits.loc[hits.index[0],"hit_def"]

        # Now go through each regular expression pattern
        e_values = []
        hit_defs = []
        paralogs = []
        for i, p in enumerate(patterns):

            for j in range(len(hits)):

                idx = hits.index[j]
                description = hits.loc[idx,"hit_def"]

                # If the pattern matches this hit description
                if p[0].search(description):

                    # If this was an NCBI blast, try to parse the NCBI line,
                    # pulling apart multi-titles
                    if ncbi_blast_db:

                        title = hits.loc[idx,"title"]
                        accession = hits.loc[idx,"accession"]

                        hd = parse_ncbi_line(title,accession)
                        if hd is not None:
                            description = hd["name"]

                    hit_defs.append(description)

                    # get e value
                    e_values.append(hits.loc[idx,"e_value"])

                    # Get paralog call
                    paralogs.append(p[1])

                    break

        # If we got at least one hit that matched
        recip_found_paralog = False
        if len(paralogs) > 0:

            # calculate posterior probability
            e_values = np.array(e_values)

            # If none of the e-values are zero, calculate the posterior
            # probability of the top hit.
            if np.sum(e_values == 0) == 0:
                pp = 1/np.array(e_values)
                pp = pp/np.sum(pp)
                match_idx = np.argsort(pp)[-1]
                match_pp = pp[match_idx]

            # If there is a zero e-value, take the first hit a zero e-value.
            # This is an (apparently) infinite posterior probaility
            else:
                match_idx = np.where(e_values == 0)[0][0]
                match_pp = np.inf

            # Get the match hit definition, paralog call, match probability
            # (relative to other possible paralogs), and difference in e value
            # relative to the best overall hit.
            recip_hit = hit_defs[match_idx]
            recip_paralog = paralogs[match_idx]
            recip_prob_match = match_pp
            recip_e_value = e_values[match_idx]

            # If the top hit has e-value of zero...
            if top_e_value == 0:

                # If the recip blast hit matches this, the difference in those
                # e-values is zero
                if recip_e_value == 0:
                    recip_del_best = 0

                # If the recip blast hit does *not* match the top, the
                # difference in those e-values is
                else:
                    recip_del_best = -np.inf
                    w = "WARNING:\n"
                    w += "BLAST returned an e-value = 0.0 for the top hit:\n\n"
                    w += f"  {hits['hit_def'].iloc[0]}\n\n"
                    w += "This hit does not match any of the patterns in paralog_patterns.\n"
                    w += "The best hit that matches something in paralog_patterns has\n"
                    w += f"an e-value of {recip_e_value} and matches '{recip_paralog}'.\n"
                    w += "These values cannot be quantitatively compared.\n"
                    w += "Topiary will set recip_found_paralog to 'False' for\n"
                    w += "this sequence.\n"
                    print(w)

            # If we get here, neither top_e_value or recip_e_value are zero.
            # (if recip_e_value was zero, it would be top_e_value and thus be
            # caught above.)
            else:
                recip_del_best = np.log10(top_e_value) - np.log10(recip_e_value)

            # Decide if we can make a paralog call based on del_best and
            # prob_match
            if recip_del_best < max_del_best and recip_prob_match > min_call_prob:
                recip_found_paralog = True

        else:
            recip_hit = top_hit_def
            recip_paralog = None
            recip_prob_match = None
            recip_del_best = None

        # Record results of analysis
        results["recip_found_paralog"].append(recip_found_paralog)
        results["recip_hit"].append(recip_hit)
        results["recip_paralog"].append(recip_paralog)
        results["recip_prob_match"].append(recip_prob_match)
        results["recip_del_best"].append(recip_del_best)

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

    num_keep_at_start = np.sum(df.keep)

    print("3",df.loc[df.loc[:,"uid"] == "tRUuOXIzEB",:])

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

    # Write out some summary statistics
    print(f"{np.sum(df.keep)} of {num_keep_at_start} sequences passed recip blast.\n")
    print("Found the following numbers of paralogs:")
    for p in patterns:
        num_of_p = np.sum(df.recip_paralog == p[1])
        print(f"    {p[1]}: {num_of_p}")
    print("",flush=True)

    return df


def recip_blast(df,
                paralog_patterns,
                local_blast_db=None,
                ncbi_blast_db=None,
                ncbi_taxid=None,
                ignorecase=True,
                max_del_best=100,
                min_call_prob=0.95,
                use_start_end=True,
                hitlist_size=50,
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
      call based on relative evalues of all paralog hits
    + `recip_del_best`: float. how much worse paralog call is than the best
      blast hit. log(e_hit_match) - log(e_hit_best)

    Parameters
    ----------

    df : pandas.DataFrame
        topiary dataframe. Will pull sequences from df.sequences. If there are
        'start' and 'stop' columns in the dataframe, only blast sequences
        between start/top (for example: start = 5, stop = 20 would blast
        sequence[5:20]. To turn this off, set use_start_end = False.
    paralog_patterns : dict
        dictionary with paralogs as values and lists of patterns to look for as
        values. Keys must be strings. Values may be strings or compiled
        re.Pattern instances. See Notes for details.
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
    max_del_best : float, default=100
        maximum value for log(e_hit_match) - log(e_hit_best) that
        allows for paralog call. This means the matched recip blast
        hit does not have to be the best hit, but must be within this
        e-value difference of the best hit. A higher number means a
        less stringent cutoff. A value of 0 would require the matched
        hit be the best hit.
    min_call_prob: float, default=0.95
        hits from all paralogs that yield a regular expression match
        are weighted by their relative e-values. Each paralog is
        assigned a relative probability. This cutoff is the minimum
        probability the best paralog match must have to result in a
        paralog call. Value should be between 0 and 1 (not inclusive),
        where min_call_prob --> 1 increases the stringency.
    use_start_end : bool, default=True
        whether or not to use start/stop columns in dataframe (if present) to
        slice subset of sequence for reciprocal blast.
    hitlist_size : int, default=100
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
    :code:`paralog_patterns` are used to match reciprocal blast hits. For
    example:

    .. code-block:: python

        paralog_patterns = {"LY96":["lymphocyte antigen 96","MD-2"],
                            "LY86":["lymphocyte antigen 86","MD-1"]}


    This would mean hits with 'lymphocyte antigen 96' or 'MD-2' will map to
    LY96; hits with 'lymphocyte antigen 86' or 'MD-1' will map to LY86. Hits
    that match neither would not be assigned a paralog.

    NOTE: string patterns are interpreted literally. The pattern :code:`"[A-Z]"`
    would be escaped to look for :code:`"\[A\-Z\]"`. If you want to use regular
    expressions in your patterns, pass them in as compiled regular expressions.
    For example,

    .. code-block:: python

        my_pattern = re.compile("[A-Z]")


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
                             max_del_best,
                             min_call_prob,
                             use_start_end)

    df, sequence_list, patterns, max_del_best, min_call_prob = out


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

    # Make paralog calls given blast output and list of patterns
    out_df = _make_recip_blast_calls(df,
                                       hit_dfs,
                                       patterns,
                                       max_del_best,
                                       min_call_prob,
                                       ncbi_blast_db)

    return out_df
