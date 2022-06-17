__description__ = \
"""
Functions for doing quality control
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-07"

from .remove_redundancy import remove_redundancy, find_cutoff
from ._clean_alignment import clean_alignment, score_alignment
from ._taxonomic_sample import taxonomic_sample
