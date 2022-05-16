__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
Reverse blast sequence datasets.
"""

from topiary import util

from .local_blast import local_blast
from .ncbi_blast import ncbi_blast
from .base import parse_ncbi_line

import pandas as pd
import numpy as np

import re, sys, copy

def _prepare_for_blast(df,
                       call_dict,
                       local_rev_blast_db,
                       ncbi_rev_blast_db,
                       ignorecase,
                       max_del_best,
                       min_call_prob,
                       use_start_end):
    """
    Check sanity of input parameters and compile patterns. Return a validated
    topiary dataframe, list of sequences, compiled set of patterns to search,
    and sanitized parameters.
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

    # Compile search patterns. Along the way, check the call_dict argument for
    # sanity
    patterns = []
    good_call_dict = True
    if type(call_dict) is not dict:
        good_call_dict = False
    else:

        # Work on copy of call dict
        call_dict = copy.deepcopy(call_dict)
        for paralog in call_dict:

            # If key isn't a string, throw error
            if type(paralog) is not str:
                good_call_dict = False
                break

            # If value isn't iterable...
            if not hasattr(call_dict[paralog],"__iter__"):

                # If it's a naked re.Pattern, put that in a length-one list
                if type(call_dict[paralog]) is re.Pattern:
                    call_dict[paralog] = [call_dict[paralog]]

                # Otherwise, throw an error
                else:
                    good_call_dict = False
                    break

            # If value is a string, put it a length-one list
            if type(call_dict[paralog]) is str:
                call_dict[paralog] = [call_dict[paralog]]

            # Deal with flags from ignorecase (should be reset fo reach paralog
            # in case user sends in compiled regex that overwrites for one
            # paralog)
            if ignorecase:
                re_flags = {"flags":re.IGNORECASE}
            else:
                re_flags = {}

            # Now go through patterns in list and prepare for regex
            to_compile = []
            for pattern in call_dict[paralog]:

                if type(pattern) is str:
                    to_compile.append(re.escape(pattern))
                elif type(pattern) is re.Pattern:
                    to_compile.append(pattern.pattern)
                    re_flags["flags"] = pattern.flags
                else:
                    good_call_dict = False
                    break

            # If we got here, looks good; compile pattern
            p = re.compile("|".join(to_compile),**re_flags)
            patterns.append((p,paralog))

    if not good_call_dict:
        err = "\ncall_dict should be a dictionary keying each paralog to a list\n"
        err += "of patterns to look for in the hit_def of each blast bit. The\n"
        err += "paralog (key) must be a key. Each pattern must be a string\n"
        err += "or a compiled re.Pattern instance. Example:\n"
        err += "    {'paralog1':['patternA','patternB',...],\n"
        err += "     'paralog2':['patternC','patternD',...],...}\n\n"
        raise ValueError(err)

    # Validate blast database arguments
    if ncbi_rev_blast_db is None and local_rev_blast_db is None:
        err = "\nPlease specificy either ncbi_rev_blast_db OR local_rev_blast_db\n\n"
        raise ValueError(err)

    if ncbi_rev_blast_db is not None and local_rev_blast_db is not None:
        err = "\nPlease specificy either ncbi_rev_blast_db OR\n"
        err += "local_rev_blast_db, but not both.\n\n"
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
    df = util.check_topiary_dataframe(df)

    # Create list of all sequences in dataframe
    sequence_list = []
    for i in range(len(df)):

        # Get sequence
        idx = df.index[i]
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
               local_rev_blast_db,
               ncbi_rev_blast_db,
               ncbi_taxid,
               hitlist_size,
               e_value_cutoff,
               gapcosts,
               local_num_threads,
               **kwargs):
    """
    Run blast on sequence_sequence list, returning a list of dataframes -- one
    df of hits for each sequence in sequence_list.
    """

    # NCBI blast
    if ncbi_rev_blast_db:

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
        w += "creating a local BLAST database for your reverse BLAST needs.\n"
        print(w)
        sys.stdout.flush()

        hit_dfs = ncbi_blast(sequence_list,
                             db=ncbi_rev_blast_db,
                             taxid=ncbi_taxid,
                             blast_program="blastp",
                             hitlist_size=hitlist_size,
                             e_value_cutoff=e_value_cutoff,
                             gapcosts=gapcosts,
                             **kwargs)

    # Local blast
    else:
        hit_dfs = local_blast(sequence_list,
                              db=local_rev_blast_db,
                              blast_program="blastp",
                              hitlist_size=hitlist_size,
                              e_value_cutoff=e_value_cutoff,
                              gapcosts=gapcosts,
                              num_threads=local_num_threads,
                              **kwargs)

    return hit_dfs


def _make_reverse_blast_calls(df,
                              hit_dfs,
                              patterns,
                              max_del_best,
                              min_call_prob,
                              ncbi_rev_blast_db):
    """
    Make paralog calls given blast output and list of patterns.
    """

    # hit_dfs is a list of dataframes, each of which has the blast hits for
    # each sequence in the input topiary dataframe.

    # Go through hits from each sequence
    results = {"reverse_found_paralog":[],
               "reverse_hit":[],
               "reverse_paralog":[],
               "reverse_prob_match":[],
               "reverse_del_best":[]}

    for blah, hits in enumerate(hit_dfs):

        # No reverse blast hits at all for this sequence
        if len(hits) == 0:
            results["reverse_found_paralog"].append(False)
            results["reverse_hit"].append(pd.NA)
            results["reverse_paralog"].append(pd.NA)
            results["reverse_prob_match"].append(np.nan)
            results["reverse_del_best"].append(np.nan)
            continue

        # Get e value of top hit
        top_e_value = hits["e_value"].iloc[0]

        # Now go through each regular expression pattern
        e_values = []
        hit_defs = []
        paralogs = []
        for i, p in enumerate(patterns):

            # Go through each hit description...
            for j, description in enumerate(hits.hit_def):

                row = hits.iloc[j]

                # If the pattern matches this hit description
                if p[0].search(description):

                    # Get hit definition
                    this_def = row["hit_def"]

                    # If this was an NCBI blast, try to parse the NCBI line,
                    # pulling apart multi-titles
                    if ncbi_rev_blast_db:
                        hd = parse_ncbi_line(row["title"],row["accession"])
                        if hd is not None:
                            this_def = hd["name"]

                    hit_defs.append(this_def)

                    # get e value
                    e_values.append(row["e_value"])

                    # Get paralog call
                    paralogs.append(p[1])

                    break


        # If we got at least one hit that matched
        reverse_found_paralog = False
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
            reverse_hit = hit_defs[match_idx]
            reverse_paralog = paralogs[match_idx]
            reverse_prob_match = match_pp
            reverse_e_value = e_values[match_idx]

            # If the top hit has e-value of zero...
            if top_e_value == 0:

                # If the reverse blast hit matches this, the difference in those
                # e-values is zero
                if reverse_e_value == 0:
                    reverse_del_best = 0

                # If the reverse blast hit does *not* match the top, the
                # difference in those e-values is
                else:
                    reverse_del_best = -np.inf
                    w = "WARNING:\n"
                    w += "\n\nBLAST returned an e-value = 0.0 for the top hit:\n\n"
                    w += f"  {hits['hit_def'].iloc[0]}\n\n"
                    w += "This hit does not match any of the patterns in call_dict.\n"
                    w += "The best hit that matches something in call_dict has\n"
                    w += f"an e-value of {reverse_e_value} and matches '{reverse_paralog}'.\n"
                    w += "Topiary will set reverse_found_paralog to 'False' for\n"
                    w += "this sequence.\n"
                    print(w)

            # If we get here, neither top_e_value or reverse_e_value are zero.
            # (if reverse_e_value was zero, it would be top_e_value and thus be
            # caught above.)
            else:
                reverse_del_best = np.log10(top_e_value) - np.log10(reverse_e_value)

            # Decide if we can make a paralog call based on del_best and
            # prob_match
            if reverse_del_best < max_del_best and reverse_prob_match > min_call_prob:
                reverse_found_paralog = True

        else:
            reverse_hit = hits["hit_def"].iloc[0]
            reverse_paralog = None
            reverse_prob_match = None
            reverse_del_best = None

        # Record results of analysis
        results["reverse_found_paralog"].append(reverse_found_paralog)
        results["reverse_hit"].append(reverse_hit)
        results["reverse_paralog"].append(reverse_paralog)
        results["reverse_prob_match"].append(reverse_prob_match)
        results["reverse_del_best"].append(reverse_del_best)

    for k in results:
        df[k] = results[k]

    # Update keep with reverse_found_paralog
    df["keep"] = np.logical_and(df["keep"],df["reverse_found_paralog"])

    return df


def reverse_blast(df,
                  call_dict,
                  local_rev_blast_db=None,
                  ncbi_rev_blast_db=None,
                  ncbi_taxid=None,
                  ignorecase=True,
                  max_del_best=100,
                  min_call_prob=0.95,
                  use_start_end=True,
                  hitlist_size=50,
                  e_value_cutoff=0.01,
                  gapcosts=(11,1),
                  local_num_threads=-1,
                  **kwargs):
    """
    Take sequences from a topiary dataframe and do a reverse blast analysis
    against an NCBI or local blast database. Looks in blast hits for the
    regular expressions defined in call_dict to call paralog for each sequence
    in df. Returns a copy of the input topiary dataframe with five new columns:

        reverse_found_paralog: True/False, whether a paralog was found
        reverse_hit: string, description for paralog hit (if found) or best hit
                     (if no match found)
        reverse_paralog: string or None. name of paralog from call_dict
        reverse_prob_match: float. probability that this is the correct paralog
                            call based on relative evalues of all paralog hits
        reverse_del_best: float. how much worse paralog call is than the best
                          blast hit. log(e_hit_match) - log(e_hit_best)

    Parameters
    ----------

    df: topiary dataframe. Will pull sequences from df.sequences. If there are
        'start' and 'stop' columns in the dataframe, only blast sequences
        between start/top (for example: start = 5, stop = 20 would blast
        sequence[5:20]. To turn this off, set use_start_end = False.

    call_dict: dictionary with paralogs as values and lists of patterns to look
               for as values. Keys must be strings. Values may be strings or
               compiled re.Pattern instances.

               example:
               {"LY96":["lymphocyte antigen 96","MD-2"],
                "LY86":["lymphocyte antigen 86","MD-1"]}

                This would mean hits with 'lymphocyte antigen 96' or 'MD-2'
                will map to LY96; hits with 'lymphocyte antigen 86' or 'MD-1'
                will map to LY86.

                NOTE: string patterns are interpreted literally. The pattern
                "[A-Z]" would be escaped to look for "\[A\-Z\]". If you want to
                use regular expressions in your patterns, pass them in as
                compiled regular expressions: pattern = re.compile("[A-Z]")

    local_rev_blast_db: local database against which to blast (incompatible with
                    ncbi_rev_blast_db)

    ncbi_rev_blast_db: database on ncbi against which to blast (incompatible
                       with local_rev_blast_db)

    ncbi_taxid: limit search to species specified by taxid for an ncbi search.
                We recommend setting this to well-annotated genomes like
                human (9606), mouse (10090), chicken (9031), and zebrafish (7955).
                This can either be a single value or a list of values.

    ignorecase: whether to ignore letter case when searching blast (default True).
                NOTE: if you pass in compiled regular expressions, the flags
                (such as re.IGNORECASE) on the *last* pattern for each paralog
                will be used for all patterns.

    max_del_best: maximum value for log(e_hit_match) - log(e_hit_best) that
                  allows for paralog call. This means the matched reverse blast
                  hit does not have to be the best hit, but must be within this
                  e-value difference of the best hit. A higher number means a
                  less stringent cutoff. A value of 0 would require the matched
                  hit be the best hit.

    min_call_prob: hits from all paralogs that yield a regular expression match
                   are weighted by their relative e-values. Each paralog is
                   assigned a relative probability. This cutoff is the minimum
                   probability the best paralog match must have to result in a
                   paralog call. Value should be between 0 and 1 (not inclusive),
                   where min_call_prob --> 1 increases the stringency.

    use_start_end: boolean. whether or not to use start/stop columns in
                    dataframe (if present) to slice subset of sequence to blast.

    histlist_size: number of hits to look at for reverse blast.

    e_value_cutoff: minimum allowable e value for a hit

    gapcosts: gap costs (must be length 2 tuple of ints)

    local_num_threads: number of threads to use for a local blast search. if -1,
                       use all available.

    kwargs: extra keyword arguments are passed directly to biopython
            NcbiblastXXXCommandline (for local blast) or qblast (for remote
            blast). These take precedence over anything specified above
            (hitlist_size, for example).
    """

    # Check sanity of input parameters and return a validated topiary dataframe,
    # list of sequences, and compiled set of patterns to search. Note: the
    # values of ncbi_rev_blast_db, local_rev_blast_db, ncbi_taxid, hitlist_size,
    # gapcosts, local_num_threads, and kwargs are checked by the blast functions
    # themselves.
    out = _prepare_for_blast(df,
                             call_dict,
                             local_rev_blast_db,
                             ncbi_rev_blast_db,
                             ignorecase,
                             max_del_best,
                             min_call_prob,
                             use_start_end)
    df, sequence_list, patterns, max_del_best, min_call_prob = out

    # Run BLAST on sequence list, returning list of dataframes -- one for each
    # seqeunce in sequence_list
    hit_dfs = _run_blast(sequence_list,
                         local_rev_blast_db,
                         ncbi_rev_blast_db,
                         ncbi_taxid,
                         hitlist_size,
                         e_value_cutoff,
                         gapcosts,
                         local_num_threads,
                         **kwargs)

    # Make paralog calls given blast output and list of patterns
    out_df = _make_reverse_blast_calls(df,
                                       hit_dfs,
                                       patterns,
                                       max_del_best,
                                       min_call_prob,
                                       ncbi_rev_blast_db)

    return out_df
