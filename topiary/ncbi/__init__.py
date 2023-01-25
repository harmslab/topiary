"""
Interface to NCBI blast and entrez.
"""

from Bio import Entrez
Entrez.email = "topiary.phylogenetics@gmail.com"

import os

# Figure out how often we can hit the NCBI servers with requests. Limit is 
# 3 times per second without an api key; 10 times per second with an api key. 
# NCBI_REQUEST_FREQ coordinates between threads so we do not exceed this limit 
try:
    Entrez.api_key = os.environ['NCBI_API_KEY']
    NCBI_REQUEST_FREQ = 1/9.5
except KeyError:
    NCBI_REQUEST_FREQ = 1/2.5

from ._parse_ncbi_line import parse_ncbi_line
from .blast import local_blast, ncbi_blast, recip_blast, make_blast_db
from .blast import records_to_df, read_blast_xml
from .blast import merge_blast_df, merge_and_annotate
from .entrez import get_sequences, get_taxid, get_proteome
