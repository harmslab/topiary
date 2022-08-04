"""
Functions for doing quality control on sequences in a dataframe.
"""

from .redundancy import remove_redundancy
from .alignment import score_alignment
from .taxonomic import get_merge_blocks

from .shrink import shrink_dataset
from .shrink import shrink_redundant
from .shrink import shrink_aligners

from .polish import polish_alignment
