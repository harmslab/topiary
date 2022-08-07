.. topiary documentation master file, created by
   sphinx-quickstart on Thu Aug 12 18:37:04 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. include:: links.rst

.. role:: red

.. role:: bolditalic

.. role:: emph

=======
topiary
=======

Python framework for doing ancestral sequence reconstruction.

+ :bolditalic:`Automatic.` Performs sequence database construction, quality
  control, multiple sequence alignment, tree construction, gene/species tree
  reconciliation, and ancestral reconstruction with minimal user input.
+ :bolditalic:`Species aware.` Integrates with the `Open Tree of Life`_
  database, improving selection of sequences and tree/ancestor inference.
+ :bolditalic:`Human-oriented.` Users prepare input as spreadsheets. Outputs are
  pdf trees, spreadsheets, clean fasta files, and graphical summaries of
  ancestor quality.
+ :bolditalic:`Flexible.` Use as a command line program or do custom analyses
  and plotting in Jupyter notebooks using the topiary API.

:bolditalic:`User input to a topiary calculation`

+------+--------------------------------------------------------------------------------------------+---------------+------------+
| name | aliases                                                                                    | species       | sequence   |
+------+--------------------------------------------------------------------------------------------+---------------+------------+
| LY96 | ESOP1;Myeloid Differentiation Protein-2;MD-2;lymphocyte antigen 96;LY-96                   | Homo sapiens  | MLPFLFF... |
+------+--------------------------------------------------------------------------------------------+---------------+------------+
| LY96 | ESOP1;Myeloid Differentiation Protein-2;MD-2;lymphocyte antigen 96;LY-96                   | Gallus gallus | MFEFVFF... |
+------+--------------------------------------------------------------------------------------------+---------------+------------+
| LY96 | ESOP1;Myeloid Differentiation Protein-2;MD-2;lymphocyte antigen 96;LY-96                   | Danio rerio   | MALWCPS... |
+------+--------------------------------------------------------------------------------------------+---------------+------------+
| LY86 | Lymphocyte Antigen 86;LY86;Myeloid Differentiation Protein-1;MD-1;RP105-associated 3;MMD-1 | Homo sapiens  | MKGFTAT... |
+------+--------------------------------------------------------------------------------------------+---------------+------------+
| LY86 | Lymphocyte Antigen 86;LY86;Myeloid Differentiation Protein-1;MD-1;RP105-associated 3;MMD-1 | Gallus gallus | MKTLNVL... |
+------+--------------------------------------------------------------------------------------------+---------------+------------+
| LY86 | Lymphocyte Antigen 86;LY86;Myeloid Differentiation Protein-1;MD-1;RP105-associated 3;MMD-1 | Danio rerio   | MKTYFNM... |
+------+--------------------------------------------------------------------------------------------+---------------+------------+

:bolditalic:`Final output tree from a topiary calculation`

.. image:: _static/img/final-tree.svg
  :align: center
  :alt: Topiary tree drawing



Quick start
===========

Installation
------------

For full installation instructions, see the

topiary can be easily installed with conda.

.. code-block:: shell-session

  conda create -n topiary topiary-asr
  conda activate topiary

Windows users and arm macOS users will need to manually install
`muscle <download_muscle_>`_ and `blast <download_blast_>`_.

Running
-------

IN PROGRESS

.. _How to cite:

How to cite
===========

If you use topiary in your research, please cite:

+ CITATION.

Please make sure to cite the tools we use in the package as well:

+ `Muscle <muscle link_>`_ (for alignment): Edgar RC (2021) *bioRxiv* https://doi.org/10.1101/2021.06.20.449169
+ `RAxML-NG`_ (for tree/ancestor inference): Kozlov et al (2019) *Bioinformatics* 35(21):4453–4455 https://doi.org/10.1093/bioinformatics/btz305
+ GeneRax_ (for gene/species tree reconcilation): Morel et al (2020) *MBE* https://doi.org/10.1093/molbev/msaa141
+ PastML_ (for ancestral gap assignment): Ishikawa et al (2019) *MBE* 36(9):2069–2085 https://doi.org/10.1093/molbev/msz131
+ `Open Tree`_ (for accessing species trees): Mctavish J, Sánchez-Reyes LL, Holster MT (2021) *Syst Biol* 70(6): 1295–1301. https://doi.org/10.1093/sysbio/syab033

Rees J & Cranston K (2017) *Biodiversity Data Journal* 5:e12581 https://doi.org/10.3897/BDJ.5.e12581

Syst Biol. 2021 Nov; 70(6): 1295–1301.
Published online 2021 May 10. doi: 10.1093/sysbio/syab033
PMCID: PMC8513759
PMID: 33970279
OpenTree: A Python Package for Accessing and Analyzing Data from the Open Tree of Life

Emily Jane Mctavish,1 Luna Luisa Sánchez-Reyes,1 and Mark T Holder2,3


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   protocol
   drawing
   data_structures
   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
