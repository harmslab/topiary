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
                         "muscle":(3,8),
                         "generax":(2,0),
                         "raxml-ng":(1,1)}


from .uid import generate_uid
from .wrap import wrap_function
from .threads import get_num_threads
