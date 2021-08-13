
__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
Core functions of topiary package, mostly for manipulation of pandas data
frames.
"""

from . import ncbi
from . import util
from . import _private

import pandas as pd
import numpy as np

# Modules for blasting, etc.
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2

from tqdm.auto import tqdm

import re, sys, os, string, random, pickle, io, urllib, http
import multiprocessing as mp

def ncbi_blast_xml_to_df(xml_files,
                         aliases=None):
    """
    xml_files: blast xml files to load. if a string, treat as a single xml file. if a
               list, treat as a list of xml files.
    aliases: dictionary for standardizing protein names.  Key specifies what
             should be output, values degenerate names that map back to that
             key.  For example:
                 "S100A9":("S100-A9","S100 A9","S-100 A9")
             would replace "S100-A9", "S100 A9", and "S-100 A9" with "S100A9"

    """

    # If only one xml file is specified, convert it to a list (of one) xml
    # file
    if type(xml_files) is str:
        xml_files = [xml_files]

    # List to hold all hits and accession numbers to download
    all_hits = []
    to_download = []

    # For each xml file
    for i, xml in enumerate(xml_files):

        # Read xml
        tmp_df = ncbi.read_blast_xml(xml)

        # Go through and
        for j in range(len(tmp_df)):

            accession = tmp_df.loc[j,"accession"]
            title = tmp_df.loc[j,"title"]
            evalue = tmp_df.loc[j,"e_value"]

            start = tmp_df.loc[j,"subject_start"]
            end = tmp_df.loc[j,"subject_end"]
            hit_info = ncbi.parse_ncbi_line(title,
                                            accession=accession,
                                            aliases=aliases)
            if hit_info is None:
                continue

            all_hits.append((xml,start,end,evalue,hit_info))
            to_download.append(hit_info["accession"])


    # Reduce to a unique set of accessions
    if len(set(to_download)) != len(to_download):
        unique_hits = []
        unique_download = []
        for i, d in enumerate(to_download):
            if d in unique_download:
                continue

            unique_hits.append(all_hits[i])
            unique_download.append(d)

        all_hits = unique_hits
        to_download = unique_download


    # Download sequences from entrez
    all_output = ncbi.entrez_download(to_download)

    # Capture sequences from the downloaded data
    captured = []
    for record in SeqIO.parse(io.StringIO(all_output), "fasta"):
        seq_id = str(record.id)
        sequence = str(record.seq)
        captured.append((seq_id,sequence))

    # Create a dictionary with appropriate keys to load into dataframe.
    out = util.create_pipeline_dict()

    # Go through every hit
    for i in range(len(all_hits)):

        # Get information from previous few rounds
        seq = captured[i][1]
        accession = captured[i][0]

        xml = all_hits[i][0]
        start = all_hits[i][1]
        end = all_hits[i][2]
        evalue = all_hits[i][3]
        hit_info = all_hits[i][4]

        # Get hit_info, if k is in key_list
        for k in hit_info.keys():
            try:
                out[k].append(hit_info[k])
            except KeyError:
                pass

        # Overwrite accession from hit_info
        out["accession"][-1] = accession

        # Load info from blast itself
        out["xml"].append(xml)
        out["sequence"].append(seq)
        out["length"].append(len(seq))
        out["evalue"].append(evalue)
        out["start"].append(start)
        out["end"].append(end)

        out["uid"].append("".join([random.choice(string.ascii_letters) for _ in range(10)]))

        out["keep"].append(True)

    return pd.DataFrame(out)

def _reverse_blast_thread(args):
    """
    Run reverse blast on a thread. Should only be called via reverse_blast.

    takes args which are interpreted as (df,i,rev_blast_db,patterns,queue)

    df: expects df has "sequence", "start", and "end" columns. Will return a
        copy of the df with "rev_hit" and "paralog" columns, corresponding
        to top hit title and call based on rev_blast_dict. It will also
        update "keep" to be False for any paralog = None sequences.
    i: iloc index to grab
    rev_blast_db: reverse blast database
    patterns: list of regular expression patterns to match sequences
    queue: multiprocessing queue for storing results
    """

    # parse args
    df = args[0]
    i = args[1]
    rev_blast_db = args[2]
    patterns = args[3]
    queue = args[4]

    index = df.index[i]

    s = df.loc[index,"sequence"]
    a = df.loc[index,"start"]
    b = df.loc[index,"end"]

    seq = s[a:b]
    hit = ncbi.local_blast(seq,db=rev_blast_db,hitlist_size=50)

    hit_call = None
    hit_def = None
    hit_e_value = None
    next_hit_e_value = None
    try:

        # Look at the top hit. See if it matches one of the patterns.
        hit_def = hit.loc[0,"hit_def"]
        for p in patterns:
            if p[0].search(hit_def):
                hit_call = p[1]
                break

        # If we made a hit call...
        if hit_call is not None:

            # Get e-value for the top hit
            hit_e_value = hit.loc[0,"e_value"]

            # Go through the next hits sequentially, looking for the first
            # hit that does not match the search pattern.  Record that
            # match e value. The ratio between hit_e_value and
            # next_hit_e_value tells us whether we are confident in the
            # reverse blast.
            for j in range(1,len(hit)):
                try:
                    next_hit_def = hit.loc[j,"hit_def"]
                    found_match = False
                    for p in patterns:
                        if p[0].search(next_hit_def):
                            found_match = True
                            break

                    # If we did not find a match, record the e-value
                    if not found_match:
                        next_hit_e_value = hit.loc[j,"e_value"]
                        break

                except KeyError:
                    # No more hits!
                    break

    except KeyError:
        pass

    # update queue
    queue.put((i,hit_def,hit_call,hit_e_value,next_hit_e_value))


def reverse_blast(df,call_dict=None,rev_blast_db="GRCh38",num_threads=-1):
    """
    df: expects df has "sequence", "start", and "end" columns. Will return a
        copy of the df with new columns:

        rev_hit: top reverse blast hit
        paralog: paralog, if pattern from call_dict matched top hit
        rev_e_value: e-value for top reverse hit
        next_rev_e_value: e-value for best reverse hit that does *not* match a
                          pattern.

        It will also update "keep" to be False for any paralog = None sequences.
    call_dict: dictionary with regular expressions as keys and calls as values.

               example:
               {"lymphogen antigen 96":"LY96",
                "MD-2":"LY96",
                "lymophogen antigen 86":"LY86",
                "MD-1":"LY86"}

    rev_blast_db: pointer to local blast database for reverse blasting.
    num_threads: number of threads to use. if -1, use all available.
    """

    print("Performing reverse blast...")

    patterns = []
    if call_dict is not None:
        for k in call_dict:
            patterns.append((re.compile(k,re.IGNORECASE),call_dict[k]))

    # Figure out number of threads to use
    if num_threads < 0:
        try:
            num_threads = mp.cpu_count()
        except NotImplementedError:
            num_threads = os.cpu_count()
            if num_threads is None:
                warning.warning("Could not determine number of cpus. Using single thread.\n")
                num_threads = 1

    # queue will hold results from each run.
    queue = mp.Manager().Queue()
    with mp.Pool(num_threads) as pool:

        # This is a bit obscure. Build a list of args to pass to the pool.
        # all_args has all len(df) reverse blast runs we want to do.
        all_args = [(df,i,rev_blast_db,patterns,queue) for i in range(len(df))]

        # Black magic. pool.imap() runs a function on elements in iterable,
        # filling threads as each job finishes. tqdm gives us a status bar.
        # By wrapping pool iterator, we get a status bar that updates as each
        # thread finishes.
        list(tqdm(pool.imap(_reverse_blast_thread,all_args),total=len(all_args)))

    # Get results out of the queue. s
    results = []
    while not queue.empty():
        results.append(queue.get())

    # Sort results
    results.sort()
    rev_hit = [r[1] for r in results]
    paralog = [r[2] for r in results]
    rev_e_value = [r[3] for r in results]
    next_rev_e_value = [r[4] for r in results]

    new_df = df.copy()
    new_df["rev_hit"] = rev_hit
    new_df["paralog"] = paralog
    new_df["rev_e_value"] = rev_e_value
    new_df["next_rev_e_value"] = next_rev_e_value

    # Remove sequences that do not reverse blast from consideration
    mask = np.array([p is None for p in paralog],dtype=np.bool)
    new_df.loc[mask,"keep"] = False

    print("Done.")

    return new_df

def remove_redundancy(df,cutoff=0.95,key_species=[]):
    """
    De-duplicate sequences according to cutoff and semi-intelligent heuristic
    criteria.

    Returns a copy of df in which "keep" is set to False for duplicates.

    This intelligently chooses between the two sequences. It favors sequences
    according to specific criteria (in this order of importance):

        1. whether sequence is from a key species
        2. whether it's annotated as structure
        3. how different the sequence length is from the median length
        4. whether it's annotated as low quality
        5. whether it's annotated as partial
        6. whether it's annotated as precursor
        7. whether it's annotated as hypothetical
        8. whether it's annotated as isoform
        9. sequence length (preferring longer)

    df: data frame with sequences.
    cutoff: %identity cutoff for combining removing sequences (0-1)
    key_species: list of key species to prefer.
    """

    def _get_quality_scores(row,key_species={}):
        """
        Get stats in order of importance (see remove_redundancy doc string).

        row: row from dataframe built by topiary
        key_species: dictionary of key species to prefer to others. only uses
                     keys for fast look up and ignores values

        returns float array for the sequence
        """

        try:
            key_species[row.species]
            key_species_score = 0.0
        except KeyError:
            key_species_score = 1.0

        return np.array([key_species_score,
                         row.structure,
                         row.diff_from_median,
                         row.low_quality,
                         row.partial,
                         row.precursor,
                         row.hypothetical,
                         row.isoform,
                         1/row.length],dtype=np.float)

    def _compare_seqs(A_seq,B_seq,A_qual,B_qual,cutoff,discard_key=False):
        """
        Compare sequence A and B based on alignment. If the sequences are
        similar within cutoff, compare A_stats and B_stats and take the sequence
        with the lower score. Scores have left-right priority.  Will select
        sequence with the first element with a lower score. If sequences are
        similar within cutoff and have equal scores, choose A.

        A_seq: sequence A
        B_seq: sequence B
        A_qual: quality scores for A
        B_qual: quality scores for B
        cutoff: cutoff for sequence comparison (~seq identity. between 0 and 1)
        discard_key: whether or not to discard key species, regardless of their
                     qualities.

        returns bool, bool

        True, True: keep both
        True, False: keep A
        False, True: keep B
        """

        # Get a normalized score: matches/len(shortest)
        score = pairwise2.align.globalxx(A_seq,B_seq,score_only=True)
        norm = score/min((len(A_seq),len(B_seq)))

        # If sequence similarity is less than the cutoff, keep both
        if norm <= cutoff:
            return True, True

        # If sequence similarity is greater than the cutoff, select one.
        else:

            # If we are not discarding key sequences and both sequences are
            # from key species, automatically keep both.
            if not discard_key:
                if A_qual[0] == 1 and B_qual[0] == 1:
                    return True, True

            # Compare two vectors. Identify first element that differs.

            # Return
            #   True, False if A is smaller
            #   False, True if B is smaller
            #   True, False if equal

            comp = np.zeros(A_qual.shape[0],dtype=np.int8)
            comp[B_qual > A_qual] = 1
            comp[B_qual < A_qual] = -1
            diffs = np.nonzero(comp)[0]

            # No difference, keep A arbitrarily
            if diffs.shape[0] == 0:
                return True, False

            # B > A at first difference, keep A
            elif comp[diffs[0]] > 0:
                return True, False

            # B < A at first difference, keep B
            else:
                return False, True


    key_species = dict([(k,None) for k in key_species])

    # This will hold output
    new_df = df.copy()

    # If not more than one seq, don't do anything
    if len(df) < 2:
        return new_df

    # Figure out how different each sequence is from the median length.  We
    # want to favor sequences that are closer to the median length than
    # otherwise.
    lengths = df.loc[df.keep,"length"]
    counts, lengths = np.histogram(lengths,bins=np.int(np.round(2*np.sqrt(len(lengths)),0)))
    median_length = lengths[np.argmax(counts)]
    new_df["diff_from_median"] = np.abs(new_df.length - median_length)

    # Get quality scores for each sequence
    quality_scores = []
    for i in range(len(new_df)):
        quality_scores.append(_get_quality_scores(new_df.iloc[i,:],key_species))

    print("Removing redundant sequences within species.")
    unique_species = np.unique(new_df.species)

    total_calcs = 0
    for s in unique_species:
        a = np.sum(new_df.species == s)
        total_calcs += a*(a - 1)//2

    with tqdm(total=total_calcs) as pbar:

        for s in unique_species:

            # species indexes are iloc row indexes corresponding to species of
            # interest.
            species_mask = new_df.species == s
            species_rows = np.arange(species_mask.shape[0],dtype=np.uint)[species_mask]
            num_this_species = np.sum(species_mask)
            for x in range(num_this_species):

                # Get species index (i). If it's already set to keep = False,
                # don't compare to other sequences.
                i = df.index[species_rows[x]]
                if not new_df.loc[i,"keep"]:
                    continue

                # Get sequence and quality score for A
                A_seq = new_df.loc[i,"sequence"]
                A_qual = quality_scores[species_rows[x]]

                # Loop over other sequences in this species
                for y in range(x+1,num_this_species):

                    # Get species index (j). If it's already set to keep = False,
                    # don't compare to other sequences.
                    j = df.index[species_rows[y]]
                    if not new_df.loc[j,"keep"]:
                        continue

                    # Get sequence and quality score for B
                    B_seq = new_df.loc[j,"sequence"]
                    B_qual = quality_scores[species_rows[y]]

                    # Decide which sequence to keep (or both). Discard sequences
                    # even if they are from key species--removing redundancy within
                    # a given species.
                    A_bool, B_bool = _compare_seqs(A_seq,B_seq,A_qual,B_qual,cutoff,
                                                   discard_key=True)

                    # Update keep for each sequence
                    new_df.loc[i,"keep"] = A_bool
                    new_df.loc[j,"keep"] = B_bool

                    # If we got rid of A, break out of this loop.  Do not need to
                    # compare to A any more.
                    if not A_bool:
                        break

            pbar.update(num_this_species*(num_this_species - 1)//2)

    print("Removing redundant sequences, all-on-all.")

    N = len(new_df)
    total_calcs = N*(N-1)//2
    with tqdm(total=total_calcs) as pbar:

        counter = 1
        for x in range(len(new_df)):

            i = new_df.index[x]

            # If we've already decided not to keep i, don't even look at it
            if not new_df.loc[i,"keep"]:
                pbar.update(N - counter)
                counter += 1
                continue

            # Get sequence of sequence and quality scores for i
            A_seq = new_df.loc[i,"sequence"]
            A_qual = quality_scores[x]

            for y in range(i+1,len(new_df)):

                j = new_df.index[y]

                # If we've already decided not to keep j, don't even look at it
                if not new_df.loc[j,"keep"]:
                    continue

                # Get sequence of sequence and quality scores for j
                B_seq = new_df.loc[j,"sequence"]
                B_qual = quality_scores[y]

                # Decide which sequence to keep (or both). Do not discard any
                # sequence from a key species.
                A_bool, B_bool = _compare_seqs(A_seq,B_seq,A_qual,B_qual,cutoff,
                                               discard_key=False)

                # Update keep for each sequence
                new_df.loc[i,"keep"] = A_bool
                new_df.loc[j,"keep"] = B_bool

                # If we got rid of A, break out of this loop.  Do not need to
                # compare to A any more.
                if not A_bool:
                    break

            pbar.update(N - counter)
            counter += 1

    print("Done.")

    return new_df

def write_fasta(df,out_file,seq_column="sequence",seq_name="pretty",
                write_only_keepers=True,empty_char="X-?",clean_sequence=False):
    """
    df: data frame to write out
    out_file: output file
    seq_column: column in data frame to use as sequence
    seq_name: column in data frame to use as >NAME.  If "pretty",
              write out a pretty names.
    write_only_keepers: whether or not to write only seq with keep = True
    empty_char: empty char. if the sequence is only empty char, do not write
                out.
    clean_sequence: replace any non-aa characters with "-"
    """

    # Make sure seq name is sane
    try:
        df[seq_name]
        take_pretty = False
    except KeyError:
        if seq_name == "pretty":
            take_pretty = True
        else:
            err = f"seq_name '{seq_name}' not recognized."
            err += "Should be a column name or 'pretty'\n"
            raise ValueError(err)

    # Make sure seq column is sane
    try:
        df[seq_column]
    except KeyError:
        err = f"seq_column '{seq_column}' not found\n."
        raise ValueError(err)

    # Construct fasta output
    out = []
    for i in range(len(df)):
        row = df.iloc[i]

        if write_only_keepers:
            if not row.keep:
                continue

        if take_pretty:
            h = _private._to_pretty(row)
        else:
            h = row[seq_name]

        seq = row[seq_column]
        is_empty = len([s for s in seq if s not in list(empty_char)]) == 0
        if seq == "" or seq is None or is_empty:
            continue

        # Replace non-aa characters with '-'
        if clean_sequence:
            seq = re.sub("[^ACDEFGHIKLMNPQRSTVWYZ-]","-",seq)

        out.append(f">{h}\n{seq}\n")

    # Write output
    f = open(out_file,"w")
    f.write("".join(out))
    f.close()


def write_phy(df,out_file,seq_column="sequence",
              write_only_keepers=True,
              empty_char="X-?",
              clean_sequence=False):
    """
    Write out a .phy file using uid as keys.

    df: data frame to write out
    out_file: output file
    seq_column: column in data frame to use as sequence
    write_only_keepers: whether or not to write only seq with keep = True
    empty_char: empty char. if the sequence is only empty char, do not write
                out.
    clean_sequence: replace any non-aa characters with "-"
    """

    # Make sure seq column is sane
    try:
        df[seq_column]
    except KeyError:
        err = f"seq_column '{seq_column}' not found\n."
        raise ValueError(err)

    if write_only_keepers:
        num_to_write = np.sum(df.keep)
    else:
        num_to_write = len(df.keep)

    all_lengths = []
    for i in range(len(df)):
        try:
            l = len(df[seq_column].iloc[i])
            if l == 0:
                raise TypeError
            all_lengths.append(l)
        except TypeError:
            if not df.keep.iloc[i] and write_only_keepers:
                pass
            else:
                err = "\n\nRow does not have sequence\n\n"
                err += f"...{df.iloc[i]}\n"
                raise ValueError(err)

    all_lengths = set(all_lengths)
    if len(all_lengths) == 0:
        err = "\n\nno sequences found\n"
        raise ValueError(err)

    if len(all_lengths) != 1:
        err = "\n\nnot all rows have the same alignment length\n\n"
        raise ValueError(err)

    # Finally, get length of alignment
    ali_length = list(all_lengths)[0]

    # Construct phy output
    out = []
    for i in range(len(df)):
        row = df.iloc[i]

        if write_only_keepers:
            if not row.keep:
                continue

        h = row["uid"]
        seq = row[seq_column]
        is_empty = len([s for s in seq if s not in list(empty_char)]) == 0
        if seq == "" or seq is None or is_empty:
            num_to_write -= 1
            continue

        # Replace non-aa characters with '-'
        if clean_sequence:
            seq = re.sub("[^ACDEFGHIKLMNPQRSTVWYZ-]","-",seq)

        out.append(f"{h}\n{seq}\n")

    out.insert(0,f"{num_to_write}  {ali_length}\n\n")

    # Write output
    f = open(out_file,"w")
    f.write("".join(out))
    f.close()

def load_fasta(df,fasta_file,load_into_column="alignment",empty_char="X-?",unkeep_missing=True):
    """
    unkeep_missing: set sequences not loading into keep = False
    empty_char: empty char. if the sequence is only empty char, set keep = False
    """

    # Create new data frame and make sure it has the column in which to load
    new_df = df.copy()
    try:
        new_df[load_into_column]
    except KeyError:
        new_df[load_into_column] = None

    # Figure out how pretty and uid calls map to df index
    pretty_to_index, uid_to_index = _private.get_index_maps(new_df)

    # Go through the fasta file and get sequences
    header = []
    seqs = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                if len(seqs) > 0:
                    seqs[-1] = "".join(seqs[-1])
                header.append(line.strip()[1:])
                seqs.append([])
            else:
                seqs[-1].extend(list(line.strip()))
    seqs[-1] = "".join(seqs[-1])

    # Load sequences from fasta into data frame
    loaded_seq = {}
    for i in range(len(header)):

        # Figure out the index to modify
        try:
            index = pretty_to_index[header[i]]
        except KeyError:
            try:
                index = uid_to_index[header[i]]
            except KeyError:
                err = f"could not map {header[i]} to data frame\n"
                raise ValueError(err)

        # Actually modify data frame
        new_df.loc[index,load_into_column] = seqs[i]

        # Record the sequence was loaded if it's not all junk (like -?X)
        if len([s for s in seqs[i] if s not in list(empty_char)]) > 0:
            loaded_seq[index] = None

    # If requested, set all sequences not in alignment to Keep = False
    if unkeep_missing:
        for i in list(new_df.index):
            try:
                loaded_seq[i]
            except KeyError:
                new_df.loc[i,"keep"] = False

    return new_df
