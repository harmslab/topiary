"""
Interface to raxml-ng.
"""

from .ancestors import generate_ancestors
from .model import find_best_model
from .tree import generate_ml_tree

from ._raxml import RAXML_BINARY
