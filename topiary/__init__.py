
# Submodules
from . import util
from . import draw
from .external import ncbi
from .external import generax
from .external import opentree
from .external import muscle

# Core functions for pipeline
from .redundancy import remove_redundancy
from .util import create_nicknames
from .external import reverse_blast
from .external import find_best_model, generate_ml_tree, generate_ancestors
from .external import reconcile

# Input/output functions
from .io import ncbi_blast_xml_to_df
from .io import write_dataframe, read_dataframe
from .io import read_fasta_into, write_fasta


def check_for_notebook():
    """
    Check whether the code is being executed in a notebook or standard
    standard python interpreter. Returns string for jupyter or IPython,
    None for something not recognized.
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
