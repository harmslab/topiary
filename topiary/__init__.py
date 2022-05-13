
from .ncbi import reverse_blast
from .io import ncbi_blast_xml_to_df, read_fasta_into, write_fasta, write_phy
from .io import write_dataframe, read_dataframe
from .redundancy import remove_redundancy
from .opentree import get_ott_id, get_species_tree
from .muscle import run_muscle
from .util import get_ott_id, create_nicknames

from . import ncbi
from . import util
from . import raxml
from . import reconcile
from . import draw

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
