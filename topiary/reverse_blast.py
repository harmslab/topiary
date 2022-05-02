__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
Reverse blast sequence datasets.
"""

from . import ncbi, util

import pandas as pd
import numpy as np

import re, sys

def reverse_blast(df,
                  call_dict,
                  local_rev_blast_db=None,
                  ncbi_rev_blast_db=None,
                  ncbi_taxid=None,
                  max_del_best=100,
                  min_call_prob=0.95,
                  use_start_stop=True,
                  hitlist_size=50,
                  e_value_cutoff=0.01,
                  gapcosts=(11,1),
                  local_num_threads=-1,
                  **kwargs):
    """
    Take sequences from a topiary dataframe and do a reverse blast analysis
    against an NCBI or local blast database. Looks in blast hits for the
    regular expressions defined in call_dict to call paralog for each sequence
    in df. Returns a new topiary dataframe with five new columns:

        reverse_found_paralog: True/False, whether a paralog was found
        reverse_hit: string, description for paralog hit (if found) or best hit
                     (if no match found)
        reverse_paralog: string or None. name of paralog from call_dict
        reverse_prob_match: float. probability that this is the correct paralog
                            call based on relative evalues of all paralog hits
        reverse_del_best: float. how much worse paralog call is than the best
                          blast hit. log(e_hit_match) - log(e_hit_best)

    df: topiary dataframe. Will pull sequences from df.sequences. If there are
        'start' and 'stop' columns in the dataframe, only blast sequences
        between start/top (for example: start = 5, stop = 20 would blast
        sequence[5:20]. To turn this off, set use_start_stop = False.

    call_dict: dictionary with paralogs as values and lists of regular
               as values. Keys and regular expressions must be strings.

               example:
               {"LY96":["lymphocyte antigen 96","MD-2"],
                "LY86":["lymphocyte antigen 86","MD-1"]}

                This would mean hits with 'lymphocyte antigen 96' and 'MD-2'
                will map to LY96; hits with 'lymphocyte antigen 86' and 'MD-1'
                will map to LY86.

    local_rev_blast_db: local database against which to blast (incompatible with
                    ncbi_rev_blast_db)

    ncbi_rev_blast_db: database on ncbi against which to blast (incompatible
                       with local_rev_blast_db)

    ncbi_taxid: limit search to species specified by taxid for an ncbi search.
                We recommend setting this to well-annotated genomes like
                human (9606), mouse (10090), chicken (9031), and zebrafish (7955).
                This can either be a single value or a list of values.

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

    use_start_stop: boolean. whether or not to use start/stop columns in
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

    # Check type of df
    if type(df) is not pd.DataFrame:
        err = "\ndf should be a topiary dataframe\n\n"
        raise ValueError(err)

    # Compile search patterns. Along the way, check the call_dict argument for
    # sanity
    patterns = []
    good_call_dict = True
    if type(call_dict) is not dict:
        good_call_dict = False
    else:
        for paralog in call_dict:

            # If key isn't a string, throw error
            if type(paralog) is not str:
                good_call_dict = False
                break

            # If value isn't iterable, thrown an error
            if not hasattr(call_dict[paralog],"__iter__"):
                good_call_dict = False
                break

            # If value is a string, make it a length-one list containing that
            # string
            if type(call_dict[paralog]) is str:
                call_dict[paralog] = [call_dict[paralog]]

            # Now go through patterns in list and make sure they are strings
            for pattern in call_dict[paralog]:
                if type(pattern) is not str:
                    good_call_dict = False
                    break

            # If we got here, looks good; compile pattern
            p = re.compile("|".join(call_dict[paralog]),re.IGNORECASE)
            patterns.append((p,paralog))

    if not good_call_dict:
        err = "\ncall_dict should be a dictionary keying each paralog to a list\n"
        err += "of regular expressions to look for in the hit_def of each blast\n"
        err += "hit. The paralog (key), as well as each regular expression,\n"
        err += "must be a string. Example: \n\n"
        err += "    {'paralog1':['regexA','regexB',...],\n"
        err += "     'paralog2':['regexC','regexD',...],...}\n\n"
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
        if max_del_best < 0:
            raise ValueError
    except (ValueError,TypeError):
        err = "\nmax_del_best must be a number between 0 and positive infinity\n\n"
        raise ValueError(err)

    # Check min_call_prob
    try:
        min_call_prob = float(min_call_prob)
        if min_call_prob <= 0 or min_call_prob >= 1:
            raise ValueError
    except (ValueError,TypeError):
        err = "\nmin_call_prob must be a number between 0 and 1 (not inclusive)\n\n"
        raise ValueError(err)

    # Check use_start_stop
    try:
        use_start_stop = bool(int(use_start_stop))
    except (ValueError,TypeError):
        err = "\nuse_start_stop must be True or False\n\n"
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
        if not use_start_stop:
            a = 0
            b = None

        # Record this sequence
        sequence_list.append(s[a:b])

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
        w = "\nBlasting against the NCBI database can be slow. You might\n"
        w += "consider creating a local BLAST database for your reverse BLAST\n"
        w += "needs.\n"
        print(w)
        sys.stdout.flush()

        hit_dfs = ncbi.ncbi_blast(sequence_list,
                                  db=ncbi_rev_blast_db,
                                  taxid=ncbi_taxid,
                                  blast_program="blastp",
                                  hitlist_size=hitlist_size,
                                  e_value_cutoff=e_value_cutoff,
                                  gapcosts=gapcosts,
                                  **kwargs)

    # Local blast
    else:
        hit_dfs = ncbi.local_blast(sequence_list,
                                   db=local_rev_blast_db,
                                   blast_program="blastp",
                                   hitlist_size=hitlist_size,
                                   e_value_cutoff=e_value_cutoff,
                                   gapcosts=gapcosts,
                                   num_threads=local_num_threads,
                                   **kwargs)

    # hit_dfs lis a list of dataframes, each of which has the blast hits for
    # each sequence in the input topiary dataframe.

    # Go through hits from each sequence
    results = {"reverse_found_paralog":[],
               "reverse_hit":[],
               "reverse_paralog":[],
               "reverse_prob_match":[],
               "reverse_del_best":[]}

    for hits in hit_dfs:

        # Get best e-value
        best_e_value = hits["e_value"].iloc[0]

        # Now go through each regular expression pattern
        e_values = []
        indexes = []
        paralogs = []
        for i, p in enumerate(patterns):

            # Go through each hit description...
            for j, description in enumerate(hits.hit_def):

                # If the pattern matches this hit description
                if p[0].search(description):

                    idx = hits.index[j]

                    # get e value, index, and paralog call for this match
                    e_values.append(hits.loc[idx,"e_value"])
                    indexes.append(idx)
                    paralogs.append(p[1])

                    break

        # If we got at least one hit that matched
        reverse_found_paralog = False
        if len(paralogs) > 0:

            # calculate posterior probability
            pp = 1/np.array(e_values)
            pp = pp/np.sum(pp)
            hit_order = np.argsort(pp)

            # Get the match hit definition, paralog call, match probability
            # (relative to other possible paralogs), and difference in e value
            # relative to the best overall hit.
            reverse_hit = hits["hit_def"].iloc[hit_order[-1]]
            reverse_paralog = paralogs[hit_order[-1]]
            reverse_prob_match = pp[hit_order[-1]]
            reverse_del_best = np.log10(best_e_value) - np.log10(e_values[hit_order[-1]])

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

    return df
