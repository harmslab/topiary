"""
Remove redundancy for datasets in a semi-intelligent way.
"""

import topiary
from topiary import check

from ._block import _check_block_redundancy

import pandas as pd
import numpy as np
from tqdm.auto import tqdm

import os
import multiprocessing as mp

# Columns to check in order
_EXPECTED_COLUMNS = ["structure","low_quality","partial","predicted",
                     "precursor","hypothetical","isoform","diff_from_median"]

class _DummyTqdm():
    """
    Fake tqdm progress bar so we don't have to show a status bar if we don't
    want to. Can be substituted wherever we would use tqdm (i.e.
    tqdm(range(10)) --> DummyTqdm(range(10)).
    """

    def __init__(self,*args,**kwargs):
        pass
    def __enter__(self):
        return self
    def __exit__(self, type, value, traceback):
        pass

    def update(self,value):
        pass

class _FakeLock():
    """
    Fake multiprocessing.Lock instance.
    """
    def __init__(self):
        pass

    def acquire(self):
        pass

    def release(self):
        pass

def _get_quality_scores(row,key_species={}):
    """
    Get stats in order of importance for a specific sequence. (see
    remove_redundancy doc string).

    Parameters
    ----------
        row: row from dataframe built by topiary
        key_species: dictionary of key species to prefer to others. only uses
                     keys for fast look up and ignores values

    Return
    ------
        float array for the sequence
    """

    try:
        if row.always_keep:
            values = [0]
        else:
            values = [1]
    except AttributeError:
        values = [1]

    # See if this is a key species
    try:
        key_species[row.species]
        values.append(0)
    except KeyError:
        values.append(1)

    # Add values for expected columns
    values.extend(list(row[_EXPECTED_COLUMNS]))

    # Flip length column
    values.append(1/len(row.sequence))

    return np.array(values,dtype=float)

def _reduce_redundancy_thread_manager(sequence_array,
                                      quality_array,
                                      keep_array,
                                      cutoff,
                                      discard_key,
                                      num_threads=-1,
                                      progress_bar=True):
    """
    Break sequence_array into a rational number of blocks given the number of
    threads and then run _check_block_redundancy on the sequence array on
    multiple threads. Updates keep_array in place.

    Parameters
    ----------
        sequence_array: array holding sequences to compare.
        quality_array: array holding vectors of quality scores, one vector for
                       each sequence
        keep_array: array holding whether or not to keep each sequence. boolean.
                    This array is updated and is the primary output of this
                    function.
        cutoff: float cutoff between 0 and 1 indicating the fractional similarity
                between two sequences above which they are considered redundant.
        discard_key: if discard_key is False, a redundant sequence will be tossed
                     even if it is from a key species
        num_threads: number of threads to use. if -1 use all available
        progress_bar: whether or not to show a progress bar

    Return
    ------
        num_threads and all_args. This is for testing/debugging purposes.
        Usually ignored. Updates keep_array in place -- this is the main output.
    """


    # Try to figure out how many cores are available
    try:
        num_machine_threads = mp.cpu_count()
    except NotImplementedError:
        num_machine_threads = os.cpu_count()

    # If we can't figure it out, revert to 1
    if num_machine_threads is None:
        print("Could not determine number of cpus. Using single thread.\n")
        num_machine_threads = 1

    # Determine number of threads useful for this problem. It's not worth
    # chopping up a super small set of comparisons
    max_useful_threads = len(sequence_array)//50
    if max_useful_threads < 1:
        max_useful_threads = 1
    if max_useful_threads > num_machine_threads:
        max_useful_threads = num_machine_threads

    # Set number of threads
    if num_threads == -1 or num_threads > max_useful_threads:
        num_threads = max_useful_threads

    # If only using one thread, don't waste overhead of making pool
    if num_threads == 1:
        block = (0,len(sequence_array))
        args = (block,block,
                sequence_array,quality_array,keep_array,
                cutoff,discard_key,_FakeLock())
        _check_block_redundancy(args)
        return 1, (args,)


    num_blocks = (num_threads**2 - num_threads)//2 + num_threads
    block_size = len(sequence_array)//num_threads

    print(f"Calculating redundancy for {num_blocks} blocks of ~{block_size} x {block_size} sequences.",flush=True)

    # Lock allowing threads to safely update keep_array
    manager = mp.Manager()
    lock =  manager.Lock()
    with mp.Pool(num_threads) as pool:

        # Calculate window sizes to cover whole L x L array, where L is number
        # of sequences
        windows = np.zeros(num_threads,dtype=int)
        windows[:] = len(sequence_array)//num_threads

        # This call spreads remainder of L/num_threads evenly across first
        # remainder windows
        windows[:(len(sequence_array) % num_threads)] += 1

        # Blocks will allow us to tile over whole redundancy matrix
        all_args = []
        for i in range(num_threads):
            for j in range(i,num_threads):

                i_block = (np.sum(windows[:i]),np.sum(windows[:i+1]))
                j_block = (np.sum(windows[:j]),np.sum(windows[:j+1]))

                all_args.append((i_block,
                                 j_block,
                                 sequence_array,
                                 quality_array,
                                 keep_array,
                                 cutoff,
                                 discard_key,
                                 lock))

        # Black magic. pool.imap() runs a function on elements in iterable,
        # filling threads as each job finishes. (Calls _check_block_redundancy
        # on every args tuple in all_args). tqdm gives us a status bar.
        # By wrapping pool.imap iterator in tqdm, we get a status bar that
        # updates as each thread finishes.
        if progress_bar:
            list(tqdm(pool.imap(_check_block_redundancy,all_args),total=len(all_args)))
        else:
            pool.imap(_check_block_redundancy,all_args)

    return num_threads, all_args


def remove_redundancy(df,
                      cutoff=0.95,
                      key_species=[],
                      silent=False,
                      only_in_species=False,
                      num_threads=-1):
    """
    Remove redundant sequences according to cutoff and semi-intelligent
    heuristics.

    Favors sequences according to the following criteria, in order:

    1. whether sequence is from a key species
    2. whether it's annotated as structure
    3. how different the sequence length is from the median length
    4. whether it's annotated as low quality
    5. whether it's annotated as partial
    6. whether it's annotated as precursor
    7. whether it's annotated as hypothetical
    8. whether it's annotated as isoform
    9. sequence length (preferring longer)

    Parameters
    ----------
    df : pandas.DataFrame
        topiary data frame with sequences
    cutoff : float, default=0.95
        %identity cutoff for combining removing sequences (between 0 and 1)
    key_species : list, default=[]
        key species to preferentially keep
    silent : bool, default=False
        whether to print output and use status bars
    only_in_species : bool, default=False
        only reduce redundancy within species; do not compare sequences between
        species
    num_threads : int, default=-1
        number of threads to use. If -1, use all available

    Returns
    -------
    topiary_dataframe : pandas.dataframe
        Copy of df in which "keep" is set to False for redundant sequences
    """

    # Process arguments
    df = check.check_topiary_dataframe(df)
    cutoff = check.check_float(cutoff,
                                           "cutoff",
                                           minimum_allowed=0,
                                           maximum_allowed=1)
    key_species = check.check_iter(key_species,
                                               "key_species",
                                               required_value_type=str,
                                               is_not_type=[str,dict])
    silent = check.check_bool(silent,"silent")

    # Encode key species as a dictionary for fast look up
    key_species = dict([(k,None) for k in key_species])

    starting_keep_number = np.sum(df.keep)

    # If not more than one seq, don't do anything
    if len(df) < 2:
        return df

    # Make sure the dataframe has the columns needed for this comparison. If
    # the dataframe does not have the column, simply set to False
    for e in _EXPECTED_COLUMNS:

        try:
            v = df[e]
            v*1.0

        # Column doesn't exist -- record as False
        except KeyError:
            df[e] = False

        # Column exists but can't be interpreted as float. Throw error.
        except TypeError:
            err = "\nThe remove_redundancy function expects a dataframe with \n"
            err += f"column '{e}' that can be interpreted as a number. This\n"
            err += f"column exists but has datatype '{df.loc[:,e].dtype}'.\n"
            err += "To fix, please rename this column and re-run this function.\n\n"
            raise ValueError(err)

    # Figure out how different each sequence is from the median length.  We
    # want to favor sequences that are closer to the median length than
    # otherwise.
    lengths = np.array([len(s) for s in df.loc[:,"sequence"]],dtype=int)
    kept_lengths = lengths[df.keep]
    hist_counts, hist_lengths = np.histogram(kept_lengths,
                                             bins=int(np.round(2*np.sqrt(len(kept_lengths)),0)))
    median_length = hist_lengths[np.argmax(hist_counts)]
    df["diff_from_median"] = np.abs(lengths - median_length)

    # Get quality scores for each sequence
    all_quality_array = []
    for i in range(len(df)):
        all_quality_array.append(_get_quality_scores(df.iloc[i,:],key_species))
    all_quality_array = np.array(all_quality_array)

    unique_species = np.unique(df.species)

    if not silent:
        P = tqdm
        print("Removing redundancy within species.",flush=True)
    else:
        P = _DummyTqdm

    with P(total=np.sum(df.keep)) as pbar:

        for i, s in enumerate(unique_species):

            species_mask = np.logical_and(df.keep,df.species == s)
            if np.sum(species_mask) < 2:
                pbar.update(np.sum(species_mask))
                continue


            # List of all sequences with keep = True as integers
            sequence_array = np.array(df.loc[species_mask,"sequence"])
            quality_array = all_quality_array[species_mask]
            keep_array = np.ones(len(sequence_array),dtype=bool)

            _reduce_redundancy_thread_manager(sequence_array=sequence_array,
                                              quality_array=quality_array,
                                              keep_array=keep_array,
                                              cutoff=cutoff,
                                              discard_key=True,
                                              num_threads=num_threads,
                                              progress_bar=False)

            # Update keep in dataframe
            df.loc[species_mask,"keep"] = keep_array

            pbar.update(np.sum(species_mask))


    if only_in_species:
        if not silent:
            final_keep_number = np.sum(df.keep)
            print(f"Reduced {starting_keep_number} --> {final_keep_number} sequences.",flush=True)
            print("Done.",flush=True)

        return df


    sequence_array = np.array(df.loc[df.keep,"sequence"])
    quality_array = all_quality_array[df.keep]
    keep_array = np.ones(len(sequence_array),dtype=bool)

    progress_bar = not silent
    if not silent:
        print("Removing redundancy in all-on-all comparison.",flush=True)

    _reduce_redundancy_thread_manager(sequence_array=sequence_array,
                                      quality_array=quality_array,
                                      keep_array=keep_array,
                                      cutoff=cutoff,
                                      discard_key=False,
                                      num_threads=num_threads,
                                      progress_bar=progress_bar)

    # Update keep array
    df.loc[df.keep,"keep"] = keep_array

    if not silent:
        final_keep_number = np.sum(df.keep)
        print(f"Reduced {starting_keep_number} --> {final_keep_number} sequences.",flush=True)
        print("Done.",flush=True)

    return df

def find_cutoff(df,
                min_cutoff=0.85,
                max_cutoff=1.00,
                try_n_values=8,
                target_number=500,
                sample_size=200,
                key_species=[]):
    """
    Find a sequence identity cutoff that yields approximately `target_number`
    sequences.

    Parameters
    ----------
    min_cutoff : float, default=0.85
        minimum identity cutoff to use (between 0 and 1)
    max_cutoff : float, default=1.00
        maximum identity cutoff to use (between 0 and 1)
    try_n_values : int, default=8
        try this many different cutoffs between min and max cutoff
    target_number : int, default=500
        find a cutoff that gets approximately this number of sequences
    sample_size : int, default=200
        grab this number of sequences from the dataframe to find the cutoff
    key_species : list, default=[]
        key species to pass to remove_redundancy

    Returns
    -------
    cutoff : float
        redundancy cutoff that yields approximately target_number sequences
    """

    df = check.check_topiary_dataframe(df)
    min_cutoff = check.check_float(min_cutoff,
                                               "min_cutoff",
                                               minimum_allowed=0,
                                               maximum_allowed=1)
    max_cutoff = check.check_float(max_cutoff,
                                               "max_cutoff",
                                               minimum_allowed=0,
                                               maximum_allowed=1)
    try_n_values = check.check_int(try_n_values,
                                               "try_n_values",
                                               minimum_allowed=2)
    target_number = check.check_int(target_number,
                                                "target_number",
                                                minimum_allowed=1)
    sample_size = check.check_int(sample_size,
                                              "sample_size",
                                              minimum_allowed=1)
    key_species = check.check_iter(key_species,
                                               "key_species",
                                               required_value_type=str,
                                               is_not_type=[str,dict])

    # Deal with cutoffs
    if min_cutoff > max_cutoff:
        err = "\nmin_cutoff must be less than or equal to max_cutoff\n\n"
        raise ValueError(err)

    # If min_cutoff and max_cutoff are the same, just return it
    if min_cutoff == max_cutoff:
        return min_cutoff
    else:
        step = (max_cutoff - min_cutoff)/(try_n_values - 1)
        cutoffs = np.array([min_cutoff + step*i for i in range(try_n_values)])

    # Make sure sample size is sane
    if sample_size > len(df.index):
        sample_size = len(df.index)

    # Create reduced dataframe to play with
    to_take = np.random.choice(df.index,sample_size,replace=False)
    small_df = df.loc[to_take,:]

    predicted_keep = []
    with tqdm(total=len(cutoffs)) as pbar:
        for c in cutoffs:
            lower_df = remove_redundancy(small_df,
                                         cutoff=c,
                                         key_species=key_species,
                                         silent=True)
            fx_kept = np.sum(lower_df.keep)/np.sum(small_df.keep)
            predicted_keep.append(fx_kept*np.sum(df.keep))

            pbar.update(1)

    predicted_keep = np.array(predicted_keep)

    # If target_number is smaller than what we got with min cutoff, return min
    # cutoff. That's as low as we can go.
    if predicted_keep[0] >= target_number:
        return min_cutoff

    # If our target number is bigger than our predicted number with the max
    # cutoff, return the maxium cutoff
    if predicted_keep[-1] <= target_number:
        return max_cutoff

    # Get the predicted_keep indexe just less than the target
    try:
        bottom = np.where(predicted_keep < target_number)[0][-1]
    except IndexError:
        bottom = None

    # Get the predicted_keep index just more than the target
    try:
        top = np.where(predicted_keep >= target_number)[0][0]
    except IndexError:
        top = None

    # Note: top and bottom should never be None at the same time because we
    # first validated that we were not entirely above the target number above.

    # If we failed to find bottom, just return top
    if bottom is None:
        return cutoffs[top]

    # If we failed to find top, just return bottom
    if top is None:
        return cutoffs[bottom]

    # Interpolate between the values just above and below the top to get cutoff
    rise = (cutoffs[top] - cutoffs[bottom])
    run = (predicted_keep[top] - predicted_keep[bottom])
    slope = rise/run
    intercept = cutoffs[top] - slope*predicted_keep[top]

    return target_number*slope + intercept
