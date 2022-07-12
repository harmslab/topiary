"""
Wrapper for generax to perform gene/species tree reconcilation.
"""

# raxml binary to use it not specified by user
GENERAX_BINARY = "generax"

import topiary
import topiary.external._interface as interface

import numpy as np

import subprocess, os, sys, time, random, string, shutil, copy
import multiprocessing as mp

def setup_generax(df,gene_tree,model,out_dir,species_tree=None):
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
    species_tree : ete3.Tree, optional
        file with species tree. if not specified, download from opentree
    """

    # -------------------------------------------------------------------------
    # Create generax data structures

    # Load gene tree
    gene_tree = topiary.io.read_tree(gene_tree)

    # Arbitrarily resolve any polytomies in gene tree
    gene_tree.resolve_polytomy()

    # Map uid to ott in the dataframe
    uid_to_ott = {}
    for i in range(len(df)):
        uid = df.uid.iloc[i]
        ott = df.ott.iloc[i]
        uid_to_ott[uid] = ott

    # For generating mapping file, record which uid are associated with
    # which ott
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


    # Make df only have uid seen (will automatically trim down to only ott
    # of interest)
    mask = np.array([u in uid_in_gene_tree for u in df.uid],dtype=np.bool)
    df = df.loc[mask]

    # Get species tree corresponding to uid seen
    if species_tree is None:
        species_tree = topiary.get_species_tree(df)

    # Resolve polytomies and make sure all branch lenghts/supports have values
    species_tree.resolve_polytomy()

    # Go across tree and set taxon names to ott, supports to 1, and branch
    # lengths to 1.
    for n in species_tree.traverse():
        if n.is_leaf():
            n.name = copy.deepcopy(n.ott)
        if n.dist != 1:
            n.dist = 1
        if n.support != 1:
            n.support = 1

    # -------------------------------------------------------------------------
    # Write out generax input

    os.mkdir(out_dir)

    # Construct the control file for generax
    control_out = []
    control_out.append("[FAMILIES]")
    control_out.append("- reconcile")
    control_out.append("starting_gene_tree = gene_tree.newick")
    control_out.append("alignment = alignment.phy")
    control_out.append("mapping = mapping.link")
    control_out.append(f"subst_model = {model}")

    # Write out control file
    f = open(os.path.join(out_dir,"control.txt"),"w")
    f.write("\n".join(control_out))
    f.close()

    # Write out .phy file
    topiary.write_phy(df,os.path.join(out_dir,"alignment.phy"),
                      seq_column="alignment")

    # Write out gene tree
    gene_tree.write(outfile=os.path.join(out_dir,"gene_tree.newick"),
                    format=5)

    # Write out species tree
    species_tree.write(outfile=os.path.join(out_dir,"species_tree.newick"),
                       format=5)

    # Write out link file
    f = open(os.path.join(out_dir,"mapping.link"),"w")
    for k in link_dict:
        f.write(f"{k}:")
        f.write(";".join(link_dict[k]))
        f.write("\n")
    f.close()


def run_generax(run_directory,
                allow_horizontal_transfer=False,
                seed=None,
                generax_binary=GENERAX_BINARY,
                log_to_stdout=True,
                other_args=[]):

    """
    Run generax. Creates a working directory, copies in the relevant files, runs
    there, and then returns to the previous directory.

    Parameters
    ----------
    run_directory : str
        directory in which to do calculation
    allow_horizontal_transfer : bool, default=False
        whether or not to allow horizontal gene transfer. This corresponds to
        the UndatedDTL (horizontal) vs UndatedDL (no horizontal) models
    seed : bool or int or str, optional
        If true, pass a randomly generated seed to generax. If int or str, use
        that as the seed (passed via --seed).
    generax_binary : str, optional
        generax binary to use
    log_to_stdout : bool, default=True
        capture log and write to std out.
    other_args : list, optional
        other arguments to pass to generax

    Returns
    -------
    generax_cmd : str
        string representation of command passed to generax
    """

    cmd = ["generax"]
    cmd.extend(["--families","control.txt"])
    cmd.extend(["--species-tree","species_tree.newick"])
    cmd.extend(["--prefix","result"])

    if allow_horizontal_transfer:
        model = cmd.extend(["--rec-model","UndatedDTL"])
    else:
        model = cmd.extend(["--rec-model","UndatedDL"])

    # seed argument is overloaded. Interpret based on type
    if seed is not None:
        if type(seed) is int:
            cmd.extend(["--seed",f"{seed:d}"])
        elif type(seed) is str:

            try:
                int(seed)
            except ValueError:
                err = f"seed {seed} could not be interpreted as an int\n"
                raise ValueError(err)

            cmd.extend(["--seed",seed])
        elif type(seed) is bool:
            if seed:
                cmd.extend(["--seed",interface.gen_seed()])
        else:
            err = "seed must be True/False, int, or string representation of int\n"
            raise ValueError(err)

    # Grab other args
    cmd.extend(other_args)

    # Make sure that generax is in the path
    try:
        subprocess.run([generax_binary])
    except FileNotFoundError:
        err = f"\ngenerax binary '{generax_binary}' not found in path\n\n"
        raise ValueError(err)

    if not os.path.exists(run_directory):
        err = f"\nrun_directory '{run_directory}' not found.\n\n"
        raise FileNotFoundError(err)

    if not os.path.isdir(run_directory):
        err = f"\nrun_directory must be a directory not found.\n\n"
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

    # Launch run
    interface.launch(cmd,run_directory,log_file)

    return " ".join(cmd)
