"""
Wrapper for generax to perform gene/species tree reconcilation.
"""

# raxml binary to use if not specified by user
GENERAX_BINARY = "generax"

import topiary
from topiary._private import interface
from topiary._private.mpi import check_mpi_configuration
from topiary._private import check

import numpy as np

import subprocess
import os
import sys
import time
import random
import string
import shutil
import copy

def _annotate_species_tree(df,species_tree,out_dir):
    """
    Clean up a species tree, making sure it is a cladogram and that the
    leaf.name attributes of the leaves are the ott. Writes to working directory.

    df : pandas.DataFrame
        topiary dataframe
    species_tree : str or ete3.Tree or dendropy.tree, optional
        species_tree to clean up and annotate. if None, pull down from the
        Open Tree of Life database.
    """

    if species_tree is not None:
        T = topiary.io.read_tree(species_tree)

    else:
        # Get species tree corresponding to uid seen
        T, dropped = topiary.df_to_species_tree(df)

        # Go across tree and set taxon names to ott, supports to 1, and branch
        # lengths to 1.
        for n in T.traverse():
            if n.is_leaf():
                n.name = copy.deepcopy(n.ott)
            if n.dist != 1:
                n.dist = 1
            if n.support != 1:
                n.support = 1

    # Resolve polytomies and make sure all branch lenghts/supports have values
    T.resolve_polytomy()

    # Make sure the leaves are formatted correctly
    for n in T.get_leaves():
        if n.name[:3] != "ott":
            err = "\nspecies tree tips must have ott labels with format ottINTEGER\n\n"
            raise ValueError(err)

    # Write out species tree
    species_tree_out = os.path.join(out_dir,"species_tree.newick")
    T.write(outfile=species_tree_out,format=5)


def _get_link_dict(df,gene_tree):
    """
    Get a dictionary linking ott (key) to a list of uid values that have that
    ott. (What genes are in a given species).

    Parameters
    ----------
    df : pandas.DataFrame
        topiary data frame
    gene_tree : ete3.Tree
        gene_tree with uid as taxon names.

    Returns
    -------
    link_dict : dict
        dictionary mapping ott to lists of uid
    uid_in_gene_tree : list
        uids seen in the gene tree
    """

    # Map uid to ott in the dataframe
    uid_to_ott = {}
    for i in df.index:
        uid = df.loc[i,"uid"]
        ott = df.loc[i,"ott"]
        uid_to_ott[uid] = ott

    # For generating mapping file, record which uid are associated with
    # which ott in gene tree we are using
    uid_in_gene_tree = []
    link_dict = {}
    for l in gene_tree.get_leaves():

        uid = l.name
        ott = uid_to_ott[uid]

        try:
            link_dict[ott].append(uid)
        except KeyError:
            link_dict[ott] = [uid]

        # Record that we saw this uid
        uid_in_gene_tree.append(uid)

    return link_dict, uid_in_gene_tree

def setup_generax(df,
                  gene_tree,
                  model,
                  out_dir,
                  keep_mask=None,
                  species_tree=None,
                  mapping_link_file=None,
                  control_file=None):
    """
    Setup a generax run directory.

    Parameters
    ----------
    df : pandas.DataFrame
        topiary data frame
    gene_tree : ete3.Tree or dendropy.Tree or str
        gene_tree with uid as taxon names. can be ete3 tree, dendropy tree, or
        newick file.
    model : str
        phylogenetic model to use (should match model used to generate gene_tree)
    out_dir : str
        output directory
    keep_mask : numpy.ndarray, optional
        bool array indicating which rows of the dataframe to keep. Note: if
        mapping_link_file is None, this mask will be ignored and recalculated.
    species_tree : str, optional
        species tree as newick file. assumes this is correctly formatted for
        a generax calculation and is copied in. if not specified, download
        species tree from opentree and annotate with ott etc
    mapping_link_file : str, optional
        mapping.link file. assumes this is correctly formatted and copied in.
        if not specified, construct.
    control_file : str, optional
        control.txt file. Assumes it is correctly formatted for a generax
        calculation and is copied in. If not specified, create.

    Returns
    -------
    keep_mask : numpy.ndarray
        bool array indicating which rows in the input datframe were kept after
        getting species etc.
    """

    # -------------------------------------------------------------------------
    # Create output directory

    os.mkdir(out_dir)

    # -------------------------------------------------------------------------
    # gene tree

    # Load gene tree, arbitrarily resolve polytomies, write out
    gene_tree = topiary.io.read_tree(gene_tree)
    gene_tree.resolve_polytomy()
    gene_tree.write(outfile=os.path.join(out_dir,"gene_tree.newick"),format=5)

    # -------------------------------------------------------------------------
    # mapping_link

    map_out = os.path.join(out_dir,"mapping.link")
    if mapping_link_file is not None:

        shutil.copy(mapping_link_file,map_out)

    else:

        link_dict, uid_in_gene_tree = _get_link_dict(df,gene_tree)

        # Write out link file
        f = open(map_out,"w")
        for k in link_dict:
            f.write(f"{k}:")
            f.write(";".join(link_dict[k]))
            f.write("\n")
        f.close()

        # Make df only have uid seen (will automatically trim down to only ott
        # of interest)
        keep_mask = df.uid.isin(uid_in_gene_tree)

    # -------------------------------------------------------------------------
    # Mask dataframe and write out alignment

    df = df.loc[keep_mask]

    # Write out .phy file
    topiary.write_phy(df,os.path.join(out_dir,"alignment.phy"),
                      seq_column="alignment",)


    # -------------------------------------------------------------------------
    # Species tree

    species_tree = _annotate_species_tree(df,species_tree,out_dir)

    # -------------------------------------------------------------------------
    # Write out generax input

    control_out_file = os.path.join(out_dir,"control.txt")
    if control_file is not None:

        shutil.copy(control_file,control_out_file)

    else:

        # Construct the control file for generax
        control_out = []
        control_out.append("[FAMILIES]")
        control_out.append("- reconcile")
        control_out.append("starting_gene_tree = gene_tree.newick")
        control_out.append("alignment = alignment.phy")
        control_out.append("mapping = mapping.link")
        control_out.append(f"subst_model = {model}")

        # Write out control file
        f = open(control_out_file,"w")
        f.write("\n".join(control_out))
        f.close()

    return keep_mask

def run_generax(run_directory,
                allow_horizontal_transfer=True,
                seed=None,
                log_to_stdout=True,
                suppress_output=False,
                other_args=[],
                write_to_script=None,
                supervisor=None,
                num_threads=1,
                generax_binary=GENERAX_BINARY):

    """
    Run generax. Creates a working directory, copies in the relevant files, runs
    there, and then returns to the previous directory.

    Parameters
    ----------
    run_directory : str
        directory in which to do calculation
    allow_horizontal_transfer : bool, default=True
        whether or not to allow horizontal gene transfer. This corresponds to
        the UndatedDTL (horizontal) vs UndatedDL (no horizontal) models
    seed : bool,int,str
        If true, pass a randomly generated seed to raxml. If int or str, use
        that as the seed. (passed via --seed)
    log_to_stdout : bool, default=True
        capture log and write to std out.
    suppress_output : bool, default=False
        whether or not to capture generax spew rather than printing to stdout.
        (ignored if log_to_stdout is True)
    other_args : list, optional
        other arguments to pass to generax
    write_to_script : str, optional
        instead of running the command, write out the command to the script file
        in the run directory. this can then be invoked later by something like
        :code:`bash script_file`.
    supervisor : Supervisor, optional
        supervisor instance to keep track of calculation inputs and outputs
    num_threads : int, default=1
        number of threads. if > 1, execute by mpirun -np num_threads
    generax_binary : str, optional
        generax binary to use

    Returns
    -------
    generax_cmd : str
        string representation of command passed to generax
    """

    if write_to_script is not None:
        if not issubclass(type(write_to_script),str):
            err = "write_to_script should be None or a string giving script name\n"
            raise ValueError(err)

    # Make sure we have full path to generax_binary
    abs_path_generax_binary = shutil.which(generax_binary)
    if abs_path_generax_binary is None:
        if supervisor is not None:
            supervisor.finalize(successful=False)
        err = f"\nrgenerax_binary '{generax_binary}' could not be found in the PATH\n\n"
        raise FileNotFoundError(err)


    if num_threads != 1:

        # Make sure we have this number of threads
        check_mpi_configuration(num_threads,abs_path_generax_binary)

        cmd = ["mpirun","-np",f"{num_threads:d}",abs_path_generax_binary]

    else:
        cmd = [abs_path_generax_binary]

    cmd.extend(["--families","control.txt"])
    cmd.extend(["--species-tree","species_tree.newick"])
    cmd.extend(["--prefix","result"])

    if allow_horizontal_transfer:
        model = cmd.extend(["--rec-model","UndatedDTL"])
    else:
        model = cmd.extend(["--rec-model","UndatedDL"])

    # seed argument is overloaded. Interpret based on type
    if seed is not None:

        # If bool and True, make the seed
        try:
            seed = check.check_bool(seed)
            if seed:
                seed = interface.gen_seed()
            else:
                seed = 0
        except ValueError:
            pass

        # Make sure the seed -- whether passed in or generated above -- is
        # actually an int.
        try:
            seed = check.check_int(seed,minimum_allowed=0)
        except ValueError:
            err = f"seed '{seed}' invalid. must be True/False or int > 0\n"
            raise ValueError(err)

        # If we have a seed > 0, append to command
        if seed > 0:
            cmd.extend(["--seed",f"{seed:d}"])

    # Grab other args
    cmd.extend(other_args)

    if not os.path.exists(run_directory):
        err = f"\nrun_directory '{run_directory}' not found.\n\n"
        raise FileNotFoundError(err)

    if not os.path.isdir(run_directory):
        err = f"\nrun_directory must be a directory.\n\n"
        raise ValueError(err)

    required_files = ["control.txt","alignment.phy",
                      "gene_tree.newick","species_tree.newick",
                      "mapping.link"]
    for f in required_files:
        filename = os.path.join(run_directory,f)
        if not os.path.exists(filename):
            err = f"\nrun_directory '{run_directory}' does not have all required\n"
            err += "files. It should have the following files:\n"
            for r in required_files:
                err += f"    {r}\n"
            err += "\n"
            raise FileNotFoundError(err)

    log_file = None
    if log_to_stdout:
        log_file = os.path.join("result","generax.log")

    if supervisor is not None:
        supervisor.event("launching generax",
                         cmd=cmd,
                         num_threads=num_threads)

    # Launch run
    try:
        interface.launch(cmd,
                         run_directory=run_directory,
                         log_file=log_file,
                         suppress_output=suppress_output,
                         write_to_script=write_to_script)
    except RuntimeError as e:
        if supervisor is not None:
            supervisor.finalize(successful=False)
        raise RuntimeError from e

    return " ".join(cmd)
