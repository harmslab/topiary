"""
Private utility functions that are not publicly exposed in the API.
"""

required_columns = ["species","name","sequence"]
reserved_columns = required_columns[:]
reserved_columns.extend(["uid","ott","alignment","keep","always_keep"])

# Data going into a newick tree can't have any of these symbols. We also reserve
# the '#' character for comments.
reserved_characters = ["(",")",";","#",":",",","'","\""]

from .uid import generate_uid
from .wrap import wrap_function
from .threads import get_num_threads
