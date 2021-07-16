
__author__ = "Michael J. Harms"
__date__ = "2021-04-08"
__description__ = \
"""
Core functions of topiary package, mostly for manipulation of pandas data
frames.
"""

from . import external
from . import util
from . import _private

import pandas as pd
import numpy as np

# Modules for blasting, etc.
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2

from tqdm.notebook import tqdm

import re, sys, os, string, random, pickle, io, urllib, http

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
        tmp_df = external.read_blast_xml(xml)

        # Go through and
        for j in range(len(tmp_df)):

            accession = tmp_df.loc[j,"accession"]
            title = tmp_df.loc[j,"title"]
            evalue = tmp_df.loc[j,"e_value"]

            start = tmp_df.loc[j,"subject_start"]
            end = tmp_df.loc[j,"subject_end"]
            hit_info = external.parse_ncbi_line(title,
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
    all_output = external.entrez_download(to_download)

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

def reverse_blast(df,call_dict=None,rev_blast_db="GRCh38"):
    """
    df: expects df has "sequence", "start", and "end" columns. Will return a
        copy of the df with "rev_hit" and "rev_call" columns, corresponding
        to top hit title and call based on rev_blast_dict. It will also
        update "keep" to be False for any rev_call = None sequences.
    call_dict: dictionary with regular expressions as keys and calls as values.

               example:
               {"lymphogen antigen 96":"LY96",
                "MD-2":"LY96",
                "lymophogen antigen 86":"LY86",
                "MD-1":"LY86"}

    rev_blast_db: pointer to local blast database for reverse blasting.
    """

    print("Performing reverse blast...")

    patterns = []
    if call_dict is not None:
        for k in call_dict:
            patterns.append((re.compile(k,re.IGNORECASE),call_dict[k]))

    rev_hit = []
    rev_call = []
    for i in tqdm(range(len(df))):

        s = df.loc[:,"sequence"].iloc[i]
        a = df.loc[:,"start"].iloc[i]
        b = df.loc[:,"end"].iloc[i]

        seq = s[a:b]
        hit = external.local_blast(seq,db=rev_blast_db,hitlist_size=1)

        hit_call = None
        try:
            hit_def = hit.loc[0,"hit_def"]
            for j in range(len(patterns)):
                if patterns[j][0].search(hit_def):
                    hit_call = patterns[j][1]
                    break
        except KeyError:
            hit_def = None

        rev_hit.append(hit_def)
        rev_call.append(hit_call)

    new_df = df.copy()
    new_df["rev_hit"] = rev_hit
    new_df["rev_call"] = rev_call

    # Remove sequences that do not reverse blast from consideration
    mask = np.array([r is None for r in rev_call],dtype=np.bool)
    new_df.loc[mask,"keep"] = False

    print("Done.")

    return new_df

def remove_redundancy(df,cutoff=0.95,key_species=[]):
    """
    De-duplicate sequences according to cutoff.

    Returns a copy of df in which "keep" is set to False for duplicates.

    This intelligently chooses between the two sequences. It favors sequences
    according to specific critera:

    If loops through these parameters, taking whichever sequence is better
    at the first parameter it encounters:

    length close to median > not low quality > not partial > not partial >
        not structure > not hypothetical > not isoform > 1/length

    """

    def _get_quality_scores(df,index,key_species={}):
        """
        Get stats in order of importance (allowed_range, low_quality, structure,
        predicted, isoform, len). These are reported such that high is
        less preferred to low.

        df: data frame with sequene data base
        index: row index (iloc)
        key_species: dictionary of key species to prefer to others. only uses
                     keys (for fast look up), ignores values

        returns new data frame with keep = False for redundant sequences removed
        from the dataset.
        """

        row = df.iloc[index]

        try:
            key_species[row.species]
            is_key_species = 0.0
        except KeyError:
            is_key_species = 1.0

        return np.array([is_key_species,
                         row.diff_from_median,
                         row.low_quality,
                         row.partial,
                         row.precursor,
                         row.structure,
                         row.hypothetical,
                         row.isoform,
                         1/row.length],dtype=np.float)

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

    print("Removing redundant sequences.")
    for i in tqdm(range(len(new_df))):

        # If we've already decided not to keep i, don't even look at it
        if not new_df.keep[i]:
            continue

        # Get sequence of sequence and quality scores for i
        A = new_df.sequence.iloc[i]
        A_stats = _get_quality_scores(new_df,i,key_species)

        for j in range(i+1,len(new_df)):

            # If we've already decided not to keep j, don't even look at it
            if not new_df.keep[j]:
                continue

            # Get sequence of sequence and quality scores for j
            B = new_df.sequence.iloc[j]
            B_stats = _get_quality_scores(new_df,i,key_species)

            # Get a normalized score: matches/len(shortest)
            score = pairwise2.align.globalxx(A,B,score_only=True)
            norm = score/min((len(A),len(B)))

            # If the sequences are similar enough
            if norm > cutoff:

                # Compare A and B stats starting from beginning and going to the end
                made_call = False
                for comp in A_stats - B_stats:

                    # If A is above B, we want B rather than A.  Set keep[i] to
                    # False and break loop.
                    if comp > 0:
                        new_df.loc[:,"keep"].iloc[i] = False
                        made_call = True
                        break

                    # If B is above A, we want A rather than B.  Set keep[j] to
                    # false and break loop.
                    elif comp < 0:
                        new_df.loc[:,"keep"].iloc[j] = False
                        made_call = True
                        break

                    # If A and B have the same score, keep comparing the next
                    # lower priority quality score
                    else:
                        continue

                # If same quality scores, toss B arbitrarily
                if not made_call:
                    new_df.loc[:,"keep"].iloc[j] = False

    print("Done.")

    return new_df

def write_fasta(df,out_file,seq_column="sequence",seq_name="pretty",
                write_only_keepers=True,empty_char="X-?"):
    """
    df: data frame to write out
    out_file: output file
    seq_column: column in data frame to use as sequence
    seq_name: column in data frame to use as >NAME.  If "pretty",
              write out a pretty names.
    write_only_keepers: whether or not to write only seq with keep = True
    empty_char: empty char. if the sequence is only empty char, do not write
                out.
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

        out.append(f">{h}\n{seq}\n")

    # Write output
    f = open(out_file,"w")
    f.write("".join(out))
    f.close()


def write_phy(df,out_file,seq_column="sequence",
              write_only_keepers=True,
              empty_char="X-?"):
    """
    Write out a .phy file using uid as keys.

    df: data frame to write out
    out_file: output file
    seq_column: column in data frame to use as sequence
    write_only_keepers: whether or not to write only seq with keep = True
    empty_char: empty char. if the sequence is only empty char, do not write
                out.
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

    ali_length = len(df[seq_column].iloc[0])

    # Construct fasta output
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
    pretty_to_index, uid_to_index = _private._get_index_maps(new_df)

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
        new_df.loc[:,load_into_column].loc[index] = seqs[i]

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
