__description__ = \
"""
Initialize topiary, exposing API.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-07"


# Submodules
from . import util
from . import draw
from . import pipeline
from . import _private
from . import _arg_processors
from . import pipeline

from .external import generax
from .external import ncbi
from .external import muscle
from .external import opentree
from .external import raxml

# Core functions for pipeline
from .quality import remove_redundancy, clean_alignment
from .util import create_nicknames
from .external import get_ott_id, get_species_tree
from .external import recip_blast
from .external import run_muscle
from .external import find_best_model, generate_ml_tree, generate_ancestors
from .external import reconcile

# Input/output functions
from .io import df_from_blast_xml, df_from_seed
from .io import write_dataframe, read_dataframe
from .io import read_fasta_into, write_fasta, write_phy

# Topiary version
from .__version__ import __version__

def check_for_notebook():
    """
    Check whether the code is being executed in a notebook or standard
    standard python interpreter.

    Return
    ------
        string for jupyter or IPython, None for something not recognized.
    """

    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return "jupyter"   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return "IPython"   # Terminal running IPython
        else:
            return None        # Not sure what interpreter

    # Probably standard Python interpreter
    except NameError:
        return None

_in_notebook = check_for_notebook()
