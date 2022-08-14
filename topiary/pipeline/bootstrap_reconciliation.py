"""
Pipeline that starts from an alignment, finds the best phylogenetic model,
builds a maximum likelihood tree, reconciles this tree with the species tree,
and then infers ancestral proteins.
"""

import topiary
from topiary.raxml import RAXML_BINARY
from topiary.generax import GENERAX_BINARY
from topiary._private import installed
from topiary._private import software_requirements
from topiary._private import check

import os
import random
import string
import shutil
import time
import json
#
# def _check_status(run_dir):
#     #XXX
#
#     complete = False
#     calc_type = None
#
#     # See if json file is there. If so, the run is done.
#     json_file = os.path.join(output,"output","run_parameters.json")
#     if os.path.isfile(json_file):
#         complete = True
#
#         f = open(json_file)
#         run_param = json.load(f)
#         f.close()
#
#         calc_type = run_param["calc_type"]
#
#     return complete, calc_type
#
# def _load_existing_pipeline(pipeline_dir):
#
#     dir_contents = [d for d in os.listdir(pipeline_dir) if os.path.isdir(d)]
#     dir_contents.sort()
#
#     # for d in dir_contents:
#     #     complete, calc_type = _check_status(d)
#     #     if complete:
#     #         if calc_type
#
#
#
#
#
#
# def alignment_to_ancestors(pipeline_dir,
#                            restart=False,
#                            overwrite=False,
#                            num_threads=-1,
#                            raxml_binary=RAXML_BINARY,
#                            generax_binary=GENERAX_BINARY):
#     """
#     Given an alignment to ancestors pipeline with ML bootstraps, calculate
#     species/gene tree reconciliation bootstraps.
#
#     Parameters
#     ----------
#     pipeline_dir : str
#         directory containing output from previous pipeline.
#     restart : bool, default=False
#         restart job from where it stopped in output directory. incompatible with
#         overwrite
#     overwrite : bool, default=False
#         whether or not to overwrite existing output. incompatible with restart
#     threads : int, default=-1
#         number of threads to use. if -1 use all available
#     raxml_binary : str, optional
#         raxml binary to use
#     generax_binary : str, optional
#         what generax binary to use
#     """
#
#     pass
#     #
#     # # --------------------------------------------------------------------------
#     # # Check sanity of overwrite, restart, and combination
#     #
#     # overwrite = check.check_bool(overwrite,"overwrite")
#     # restart = check.check_bool(restart,"restart")
#     #
#     # if overwrite and restart:
#     #     err = "overwrite and restart flags are incompatible.\n"
#     #     raise ValueError(err)
#     #
#     # num_threads = check.check_int(num_threads,
#     #                               "num_threads",
#     #                               minimum_allowed=-1)
#     # if num_threads == 0:
#     #     err = "num_threads should be -1 (use all available) or a integer > 0\n"
#     #     err += "indicating the number of threads to use.\n"
#     #     raise ValueError(err)
#     #
#     # # --------------------------------------------------------------------------
#     # # Validate software stack required for this pipeline
#     #
#     # to_validate = [{"program":"raxml-ng",
#     #                 "binary":raxml_binary,
#     #                 "min_version":software_requirements["raxml-ng"],
#     #                 "must_pass":True}]
#     # to_validate.append({"program":"generax",
#     #                     "binary":generax_binary,
#     #                     "min_version":software_requirements["generax"],
#     #                     "must_pass":True})
#     # to_validate.append({"program":"mpirun",
#     #                     "min_version":software_requirements["mpirun"],
#     #                     "must_pass":True})
#     #
#     # installed.validate_stack(to_validate)
#     #
#     # counter = 0 # HACK
#     #
#     #
#     # output = f"{counter:02d}_reconciliation-bootstraps"
#     # run_calc = _check_restart(output,restart)
#     # if run_calc:
#     #
#     #     topiary.reconcile(previous_dir=previous_dir,
#     #                       output=output,
#     #                       allow_horizontal_transfer=allow_horizontal_transfer,
#     #                       generax_binary=generax_binary,
#     #                       num_threads=num_threads,
#     #                       bootstrap=True)
#     # counter += 1
#     #
#     # # Generate final ancestors on tree with branch supports
#     # previous_dir = output
#     # output = f"{counter:02d}_ancestors_with_branch_supports"
#     # run_calc = _check_restart(output,restart)
#     # if run_calc:
#     #     topiary.generate_ancestors(previous_dir=previous_dir,
#     #                                output=output,
#     #                                num_threads=1,
#     #                                alt_cutoff=alt_cutoff)
#     #
#     # os.chdir(current_dir)
