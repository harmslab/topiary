"""
Private utility functions that are not publicly exposed in the API.
"""

required_columns = ["species","name","sequence"]
reserved_columns = required_columns[:]
reserved_columns.extend(["uid","ott","alignment","keep","always_keep"])

# Data going into a newick tree can't have any of these symbols. We also reserve
# the '#' character for comments.
reserved_characters = ["(",")",";","#",":",",","'","\""]

# External (non-python) software version requirements
software_requirements = {"blastp":(2,0),
                         "makeblastdb":(2,0),
                         "muscle":(5,0),
                         "generax":(2,0),
                         "raxml-ng":(1,1),
                         "mpirun":(0,0)}


from .uid import generate_uid
from .supervisor import Supervisor
