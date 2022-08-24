"""
Final step on the pipeline. Run replicates in embarassingly parallel fashion
across compute nodes using MPI.
"""

import topiary
from topiary.raxml import RAXML_BINARY
from topiary.generax import GENERAX_BINARY
from topiary.generax._reconcile_bootstrap import reconcile_bootstrap
from topiary._private import installed
from topiary._private import software_requirements
from topiary._private.mpi import check_mpi_configuration
from topiary._private import check
from topiary._private import Supervisor

import os
import datetime
import glob
import shutil


def bootstrap_reconcile(previous_run_dir,
                        num_threads,
                        restart=False,
                        overwrite=False,
                        raxml_binary=RAXML_BINARY,
                        generax_binary=GENERAX_BINARY):
    """
    Perform a bootstrap branch support calculation for the gene/species tree
    reconciliation portion of the analysis.

    previous_run_dir : str
        previous pipeline run directory. should end with last directory added
        as 04_bootstraps.
    num_threads : int
        number of threads to use. GeneRax uses MPI for parallelization;
        num_threads correspond to the number of "slots" in MPI lingo. This job
        can be massively parallelized as there is no cross-talk between
        replicates, so feel free to span this across as many compute nodes as
        you like.
    restart : bool, default=False
        restart job from where it stopped in output directory. incompatible with
        overwrite
    overwrite : bool, default=False
        whether or not to overwrite existing output. incompatible with restart.
        This will overwrite an existing 05_reconciliation-bootstraps directory,
        not the rest of the pipeline directory.
    threads : int, default=-1
        number of threads to use. if -1 use all available
    raxml_binary : str, optional
        raxml binary to use
    generax_binary : str, optional
        what generax binary to use
    """

    # Make sure pipeline directory is present
    if not os.path.isdir(previous_run_dir):
        err = f"\nprevious_run_dir '{previous_run_dir}' does not exist\n\n"
        raise FileNotFoundError(err)

    # --------------------------------------------------------------------------
    # Check sanity of num_threads

    num_threads = check.check_int(num_threads,
                                  "num_threads",
                                  minimum_allowed=1)

    # --------------------------------------------------------------------------
    # Check sanity of overwrite, restart, and combination

    overwrite = check.check_bool(overwrite,"overwrite")
    restart = check.check_bool(restart,"restart")

    if overwrite and restart:
        err = "overwrite and restart flags are incompatible.\n"
        raise ValueError(err)

    # --------------------------------------------------------------------------
    # Validate software stack required for this pipeline

    to_validate = [{"program":"raxml-ng",
                    "binary":raxml_binary,
                    "min_version":software_requirements["raxml-ng"],
                    "must_pass":True}]

    to_validate.append({"program":"generax",
                        "binary":generax_binary,
                        "min_version":software_requirements["generax"],
                        "must_pass":True})

    to_validate.append({"program":"mpirun",
                        "min_version":software_requirements["mpirun"],
                        "must_pass":True})

    installed.validate_stack(to_validate)

    # If we got here, reconciliation software is ready to go. Now check to
    # whether mpi can really grab the number of slots requested.
    check_mpi_configuration(num_threads,generax_binary)

    # --------------------------------------------------------------------------
    # Validate the previous calculation

    os.chdir(previous_run_dir)

    if not os.path.isdir("04_bootstraps"):
        err = f"previous_run_dir '{previous_run_dir}' does not have an 04_bootstraps\n"
        err += "directory. This directory is necessary as the input to the a\n"
        err += "reconciliation bootstrap calculation.\n\n"
        os.chdir("..")
        raise FileNotFoundError(err)

    # Load calculation and make sure it completed
    supervisor = Supervisor("04_bootstraps")
    if supervisor.status != "complete":
        err = f"{previous_run_dir}/04_bootstraps exists but has status '{supervisor.calc_status}'\n"
        if supervisor.status == "empty":
            err += "It does not appear this calculation has been run.\n\n"
        elif supervisor.status == "running":
            err += "This job is either still running or crashed.\n\n"
        else:
            err += "This job crashed before completing.\n\n"
        os.chdir("..")
        raise RuntimeError(err)

    # Get number of replicates. Make sure user did not request more slots than
    # replicates.
    num_replicates = len(glob.glob(os.path.join("04_bootstraps",
                                                "output",
                                                "bootstrap_replicates",
                                                "*.phy")))
    if num_threads > num_replicates:
        err = f"The number of requested threads (slots: {num_threads}) is\n"
        err += f"greater than the number of bootstrap replicates ({num_replicates}).\n"
        err += f"Please restart the calculation requesting no more than {num_replicates} slots.\n\n"
        os.chdir("..")
        raise ValueError(err)


    # Try to get a time estimate from supervisor stack
    try:
        for p in supervisor.previous_entries:
            if p["calc_type"] == "reconcile":

                prev_dt = p["completion_time"] - p["creation_time"]
                prev_num_threads = p["events"][1]["num_threads"]

                pretty_prev = f"{str(datetime.timedelta(seconds=prev_dt))} (D:H:M:S)"

                t_per_rep = prev_dt*prev_num_threads
                t_overall = t_per_rep*num_replicates
                t_per_slot = t_overall/num_threads

                pt = f"{str(datetime.timedelta(seconds=t_per_slot))} (D:H:M:S)"

                out = ["\n----------------------------------------------------------------------\n"]
                out.append(f"The first reconciliation calculation took {pretty_prev}")
                out.append(f"to complete on {prev_num_threads} threads. Assuming a similar machine architecture,")
                out.append("topiary predicts the current calculation should take ")
                out.append(f"{pt} to complete on {num_threads} threads.")
                out.append("\n----------------------------------------------------------------------\n")
                print("\n".join(out))

    except (KeyError,IndexError):
        out = ["\n----------------------------------------------------------------------\n"]
        out.append("Could not determine previous run time and thread/counts to")
        out.append("estimate run time for this calculation.")
        out.append("\n----------------------------------------------------------------------\n")
        print("\n".join(out))

    # Make sure the output either exists with --overwrite or --restart or
    # does not exist
    calc_dir = "05_reconcile-bootstraps"
    if os.path.isdir(calc_dir):
        if overwrite:
            shutil.rmtree(calc_dir)

        if not restart:
            err = f"'{previous_run_dir}/{calc_dir}' already exists. Either remove\n"
            err += "it, specify --overwrite, or specify --restart.\n\n"
            os.chdir("..")
            raise FileExistsError(err)

    if os.path.isdir(calc_dir) and restart:

        if not os.path.isdir(os.path.join(calc_dir,"working","replicates")):
            err = "\ncould not restart the calculation. Please delete the directory\n"
            err += f"'{previous_run_dir}/{calc_dir}' and try again.\n\n"
            os.chdir("..")
            raise ValueError(err)

        supervisor = Supervisor(calc_dir)
        supervisor.update("calc_status","running")
        supervisor.event("Restarting calculation.")

        reconcile_bootstrap(supervisor.df,
                            supervisor.model,
                            supervisor.gene_tree,
                            supervisor.species_tree,
                            supervisor.reconciled_tree,
                            supervisor.run_parameters["allow_horizontal_transfer"],
                            supervisor.seed,
                            bootstrap_directory=None,
                            restart="replicates",
                            overwrite=False,
                            supervisor=supervisor,
                            num_threads=num_threads,
                            generax_binary=generax_binary,
                            raxml_binary=raxml_binary)

    else:

        topiary.reconcile(prev_calculation="04_bootstraps",
                          calc_dir=calc_dir,
                          bootstrap=True,
                          overwrite=False,
                          num_threads=num_threads,
                          raxml_binary=raxml_binary,
                          generax_binary=generax_binary)

    os.chdir('..')
