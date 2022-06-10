__description__ = \
"""
Interface to NCBI blast and entrez.
"""
__author__ = "Michael J. Harms"
__date__ = "2022-06-07"

from Bio import Entrez
Entrez.email = "topiary.phylogenetics@gmail.com"

from . _parse_ncbi_line import parse_ncbi_line
from . blast import read_blast_xml, local_blast, merge_blast_df, ncbi_blast, recip_blast, make_blast_db
from . entrez import get_sequences, get_taxid, get_proteome
