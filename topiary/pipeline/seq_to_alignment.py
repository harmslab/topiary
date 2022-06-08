
import topiary
from topiary import _arg_processors
import pandas as pd
import numpy as np

import os, re, shutil

from .base import _run_and_print


def read_paralog_patterns(paralog_patterns):

    err_summary = ["This file should have the following format:\n",
                   "paralog1;pattern1;pattern2;pattern3"
                   "paralog2;pattern4;pattern5;pattern6",
                   "...\n",
                   "For a real pair of proteins, this could be:\n",
                   "LY96;MD2;MD-2;ESOP1",
                   "LY86;MD1;MD-1\n",
                   "Blank lines and anything after # on a line are ignored\n"]
    err_summary = "\n".join(err_summary)

    if type(paralog_patterns) is str:
        if os.path.isfile(paralog_patterns):

            out_dict = {}
            with open(paralog_patterns) as f:
                for line in f:

                    # Remove everything after '#' (comment)
                    l = line.split("#")[0].strip()

                    # Skip blank lines
                    if l.strip() == "":
                        continue

                    # Now split line on ";" and get non-"" entries.
                    col = [c.strip() for c in l.split(";")]
                    col = [c for c in col if c != ""]

                    if len(col) < 2:
                        err = f"\nCould not parse line '{line}' in paralog_patterns.\n\n"
                        err += err_summary
                        raise ValueError(err)

                    # Create pattern dictionary, searching for all unique patterns
                    out_dict[col[0]] = list(set(col))

            paralog_patterns = out_dict

    patterns = _arg_processors.process_paralog_patterns(paralog_patterns)

    return dict([(p[1],p[0]) for p in patterns])


def read_key_species(key_species):
    """
    """

    err_summary = ["This file should have the following format:\n",
                   "species 1;species 2;species 3; ...\n\n"
                   "For example, this could be:\n",
                   "Homo sapiens;Mus musculus;Danio rerio\n\n",
                   "Blank lines and anything after # on a line are ignored\n"]
    err_summary = "\n".join(err_summary)

    if type(key_species) is str:
        if os.path.isfile(key_species):

            out_list = []
            with open(key_species) as f:
                for line in f:

                    # Remove everything after '#' (comment)
                    l = line.split("#")[0].strip()

                    # Skip blank lines
                    if l.strip() == "":
                        continue

                    # Now split line on ";" and get non-"" entries.
                    col = [c.strip() for c in l.split(";")]
                    col = [c for c in col if c != ""]

                    out_list.extend(col)

            if len(out_list) == 0:
                w = f"\nWarning: no key species found in '{key_species}'\n"
                print(w,flush=True)

            key_species = out_list[:]

    key_species = _arg_processors.process_iter(key_species,
                                               "key_species",
                                               required_value_type=str,
                                               is_not_type=[str,dict])

    return key_species

def read_ncbi_taxid(ncbi_rev_blast_taxid):

    taxid_out = []
    taxid_is_bad = True

    # Is it int-like? Convert integer input to list of one string. We're
    # stringent on this -- no floats -- because an int cast will silently
    # round down. If someone used a taxid like 960.6 (extra .) it would be
    # interpreted as 960, which is very much the wrong taxid
    if np.issubdtype(type(ncbi_rev_blast_taxid),np.integer):
        ncbi_rev_blast_taxid = [f"{ncbi_rev_blast_taxid}"]

    # This wacky line sees if the taxid is iterable, but not a type *class*.
    # Catches weird edge case where user passes in str or int as a taxid
    if hasattr(ncbi_rev_blast_taxid,"__iter__") and type(ncbi_rev_blast_taxid) is not type:

        taxid_is_bad = False

        # If taxid is a single string, convert it to a list of one string. Split
        # on "," in case coming in from command line.
        if type(ncbi_rev_blast_taxid) is str:
            ncbi_rev_blast_taxid = [t.strip() for t in ncbi_rev_blast_taxid.split(",")]

        # Go through list of taxids and put in correct format
        for t in ncbi_rev_blast_taxid:
            if type(t) is str:
                taxid_out.append(f"txid{t}[ORGN]")
            elif np.issubdtype(type(t),np.integer):
                taxid_out.append(f"txid{t}[ORGN]")
            else:
                taxid_is_bad = True
                break

    # If taxid was None to begin with, ignore it
    else:
        if ncbi_rev_blast_taxid is None:
            taxid_out = []
            taxid_is_bad = False

    if taxid_is_bad:
        err = "\ntaxid should be either a single ncbi taxon id (e.g. 9606 for\n"
        err += "human) or list of ncbi taxon ids. Individual ids can either be\n"
        err += "int or str. These are passed to NCBI without further validation.\n\n"
        raise ValueError(err)


def rockit(xml_input,
           paralog_patterns,
           key_species=[],
           output_dir="rockit",
           min_redundancy_cutoff=0.85,
           target_number_of_sequences=500,
           ncbi_rev_blast_db="nr",
           local_rev_blast_db=None,
           ncbi_rev_blast_taxid="9606",
           phylo_context="All life",
           clean_alignment=True,
           overwrite=False):
    """

    Parameters
    ----------
        xml_input: xml file, directory containing xml files, or list of xml files

        paralog_patterns: dictionary of paralog patterns or file containing
                          paralog_patterns

        key_species: list of key species or file containing key species

        output_dir: output directory

        min_redundancy_cutoff: when selecting a redundancy cutoff, do not select
                               a cutoff less than this

        target_number_of_sequences: find a redundancy cutoff that yields
                                    approximately this number of sequences in
                                    final alignment

        ncbi_rev_blast_db: NCBI reverse blast database. Ignored if local_rev_blast_db
                           is specified

        local_rev_blast_db: local reverse blast database to use.

        ncbi_rev_blast_taxid: limit NCBI reverse blast to the specified taxid.
                              integer or list of integers. If called from the
                              command line use the format:
                              --ncbi_rev_blast_db=9606,10090

        phylo_context: look for species within a specific phylogenetic context
                       on open tree of life. To see what's available given the
                       current opentree database, type the following on your
                       python terminal:

                            from opentree import OT
                            print(OT.tnrs_contexts().response_dict)

        clean_alignment: whether or not to try to clean up the alignment to
                         remove sparse/gappy sequences.

        overwrite: overwrite outputs. True or False. On command line, use
                    --overwrite.

    Return
    ------
        None.

        Creates and populates output_dir with output files. Key outputs are
        final-dataframe.csv and final-alignment.fasta. All files start with
        pipeline step number.
    """

    # Get current workign directory
    current_dir = os.getcwd()

    # See if this is an input directory. If so, make it a full absolute path
    if os.path.exists(str(xml_input)):
        xml_input = os.path.abspath(str(xml_input))

    # Deal with the paralog_patterns argument (either read file or just load dict)
    paralog_patterns = read_paralog_patterns(paralog_patterns)

    # Deal with the key_species argument (either read file or just load dict)
    key_species = read_key_species(key_species)

    # check min_redundancy_cutoff
    min_redundancy_cutoff = _arg_processors.process_float(min_redundancy_cutoff,
                                                          "min_redundancy_cutoff",
                                                          minimum_allowed=0,
                                                          maximum_allowed=1)

    # check target_number_of_sequences
    target_number_of_sequences = _arg_processors.process_int(target_number_of_sequences,
                                                             "target_number_of_sequences",
                                                             minimum_allowed=1)

    # If local_rev_blast_db is not None, create absolute path to that file
    if local_rev_blast_db is not None:
        local_rev_blast_db = str(local_rev_blast_db)
        local_rev_blast_db = os.path.join(current_dir,local_rev_blast_db)

        # Check for existance of blast database
        blast_file = f"{local_rev_blast_db}.psq"
        if not os.path.isfile(blast_file):
            err = f"\ncould not find local blast database {local_rev_blast_db}\n\n"
            raise ValueError(err)

        # Set ncbi reverse blast to None -- override with local
        ncbi_rev_blast_db = None

    # Make ncbi_rev_blast_db a string
    if ncbi_rev_blast_db is not None:
        ncbi_rev_blast_db = str(ncbi_rev_blast_db)
        ncbi_rev_blast_taxid = read_ncbi_taxid(ncbi_rev_blast_taxid)

    # Make sure the phylogenetic context is sane
    phylo_context = topiary.opentree.is_allowed_phylo_context(phylo_context)

    # Check for existance of output_directory
    output_dir = str(output_dir)
    if os.path.isdir(output_dir):
        if overwrite:
            shutil.rmtree(output_dir)
        else:
            err = f"\noutput directory '{output_dir}' already exists.\n"
            err += "To overwrite, set overwrite = True.\n\n"
            raise FileExistsError(err)

    # Move into output directory.
    current_dir = os.getcwd()
    os.mkdir(output_dir)
    os.chdir(output_dir)

    step_counter = 0

    # Load xml files and download sequences
    df, step_counter = _run_and_print(function=topiary.ncbi_blast_xml_to_df,
                                      kwargs={"xml_input":xml_input},
                                      step_counter=step_counter,
                                      out_file_string="initial",
                                      human_string="Loading xml files.")

    # Make human readable sequence nicknames
    df, step_counter = _run_and_print(function=topiary.create_nicknames,
                                      kwargs={"df":df,"paralog_patterns":paralog_patterns},
                                      step_counter=step_counter,
                                      out_file_string="nickname",
                                      human_string="Assigning sequence nicknames.")

    # Get ott id. We do this early on in case we can't find a species ott. This
    # will set those sequences to False
    df, step_counter = _run_and_print(function=topiary.get_ott_id,
                                      kwargs={"df":df,"phylo_context":phylo_context},
                                      step_counter=step_counter,
                                      out_file_string="ott",
                                      human_string="Getting OTT (species) ids.")


    # Remove basically identical sequences early on. Goal is to minimize
    # sequences passed to reverse blast.
    df, step_counter = _run_and_print(function=topiary.remove_redundancy,
                                      kwargs={"df":df,"cutoff":0.999,"key_species":key_species,"only_in_species":True,"silent":False},
                                      step_counter=step_counter,
                                      out_file_string="remove-identical",
                                      human_string="Removing identical sequences within species.")


    # If a reverse blast database was passed in, do reverse blast.
    if local_rev_blast_db is not None or ncbi_rev_blast_db is not None:

        kwargs = {"df":df,
                  "paralog_patterns":paralog_patterns,
                  "local_rev_blast_db":local_rev_blast_db,
                  "ncbi_rev_blast_db":ncbi_rev_blast_db,
                  "ncbi_taxid":ncbi_rev_blast_taxid}
        df, step_counter = _run_and_print(function=topiary.reverse_blast,
                                          kwargs=kwargs,
                                          step_counter=step_counter,
                                          out_file_string="reverse-blast",
                                          human_string="Running reverse blast.")

    # Identify a redundancy cutoff that yields approximately the target number
    # of sequences.
    print("-------------------------------------------------------------------")
    print(f"Identifying redundancy cutoff that yields ~ {target_number_of_sequences} sequences.")
    print("-------------------------------------------------------------------")
    print("",flush=True)
    cutoff = topiary.quality.find_cutoff(df,
                                         min_cutoff=min_redundancy_cutoff,
                                         target_number=target_number_of_sequences)
    print()
    print(f"Found cutoff: {cutoff:.2f}\n",flush=True)

    # Lower sequence redundancy
    df, step_counter = _run_and_print(function=topiary.remove_redundancy,
                                      kwargs={"df":df,"cutoff":cutoff,"key_species":key_species},
                                      step_counter=step_counter,
                                      out_file_string="lower-redundancy",
                                      human_string="Lowering sequence redundancy.")



    # Clean alignment
    if clean_alignment:

        df, step_counter = _run_and_print(function=topiary.run_muscle,
                                          kwargs={"input":df},
                                          step_counter=step_counter,
                                          out_file_string="first-pass-alignment",
                                          human_string="Doing initial alignment.")

        df, step_counter = _run_and_print(function=topiary.quality.clean_alignment,
                                          kwargs={"df":df,"key_species":key_species},
                                          step_counter=step_counter,
                                          out_file_string="cleaned-up-alignment",
                                          human_string="Removing difficult-to-align sequences.")

    # Align final set of high quality sequences
    df, step_counter = _run_and_print(function=topiary.run_muscle,
                                      kwargs={"input":df},
                                      step_counter=step_counter,
                                      out_file_string="final-alignment",
                                      human_string="Doing final alignment.")

    alignment_file = f"{step_counter:02d}_final-alignment.fasta"

    print("-------------------------------------------------------------------")
    print(f"Writing final alignment to {os.path.join(output_dir,alignment_file)}")
    print("-------------------------------------------------------------------")
    print("",flush=True)

    topiary.write_fasta(df,alignment_file,
                        seq_column="alignment")

    os.chdir(current_dir)
