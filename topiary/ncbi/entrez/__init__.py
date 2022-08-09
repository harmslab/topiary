"""
Functions for interfacing with NCBI entrez databases.
"""

from .proteome import get_proteome
from .proteome import get_proteome_ids
from .download import ncbi_ftp_download
from .taxid import get_taxid
from .sequences import get_sequences
