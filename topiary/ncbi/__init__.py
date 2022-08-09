"""
Interface to NCBI blast and entrez.
"""

from Bio import Entrez
Entrez.email = "topiary.phylogenetics@gmail.com"

from ._parse_ncbi_line import parse_ncbi_line
from .blast import local_blast, ncbi_blast, recip_blast, make_blast_db
from .blast import records_to_df, read_blast_xml
from .blast import merge_blast_df, merge_and_annotate
from .entrez import get_sequences, get_taxid, get_proteome
