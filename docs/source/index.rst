.. topiary documentation master file, created by
   sphinx-quickstart on Thu Aug 12 18:37:04 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. role:: red

.. role:: bolditalic

.. role:: emph

=======
topiary
=======

Python framework for doing ancestral sequence reconstruction.

+ Automates sequence database construction, quality control, multiple
  sequence alignment, tree construction, gene/species tree reconciliation,
  and ancestral reconstruction.
+ Spreadsheet centric. Users prepare input as spreadsheets. The sequence database
  and alignments are stored as csv files.
+ Final outputs are human-readable trees and ancestral sequences.
+ Use with simple command line tools or do custom analyses using the topiary
  API.


.. _Installation:

Installation
============

------------
Requirements
------------

+ Core scientific python libraries:

  + `Python >= 3.8 <https://www.python.org/downloads/>`_
  + `numpy <https://numpy.org>`_
  + pandas_

+ Tree manipulation/drawing:

  + `ete3 <http://etetoolkit.org/download/>`_
  + `toytree <https://toytree.readthedocs.io/en/latest/3-installation.html>`_

+ Packages used for tree/ancestor inferences:

  + `NCBI BLAST+ >= 2.2 <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_
  + `muscle >= 5.0 <https://github.com/rcedgar/muscle/releases/>`_
  + `GeneRax >= 2.0 <https://github.com/BenoitMorel/GeneRax/releases/>`_
  + `RAxML-NG >= 1.1 <https://github.com/amkozlov/raxml-ng/releases/>`_
  + `pastml >= 1.9 <https://pastml.pasteur.fr>`_


---------------
macOS and linux
---------------

We recommend installing topiary via miniconda_.

.. code-block:: shell-session

  conda install -c bioconda topiary-asr


On a command line/conda terminal:

.. code-block:: shell-session

  conda install -c bioconda muscle
  conda install -c bioconda blast


If you plan to run generax and raxml-ng, and are on macOS or linux:

.. code-block:: shell-session

  conda install -c bioconda generax
  conda install -c bioconda raxml-ng


Finally:

.. code-block:: shell-session

  git clone https://github.com/harmslab/topiary.git
  cd topiary
  python setup.py install


To make sure the code installed properly, run:

.. code-block:: shell-session

  conda install pytest
  pytest


.. _How to cite:

How to cite
===========

If you use topiary in your research, please cite:
+ CITATION.

Please make sure to cite the tools we use in the package as well:

+ Muscle_ (for alignment): Edgar RC (2021) *bioRxiv* https://doi.org/10.1101/2021.06.20.449169
+ `RAxML-NG`_ (for tree/ancestor inference): Kozlov et al (2019) *Bioinformatics* 35(21):4453–4455 https://doi.org/10.1093/bioinformatics/btz305
+ GeneRax_ (for gene/species tree reconcilation): Morel et al (2020) *MBE* https://doi.org/10.1093/molbev/msaa141
+ PastML_ (for ancestral gap assignment): Ishikawa et al (2019) *MBE* 36(9):2069–2085 https://doi.org/10.1093/molbev/msz131
+ `Open Tree of Life`_ (for species trees): Rees J & Cranston K (2017) *Biodiversity Data Journal* 5:e12581 https://doi.org/10.3897/BDJ.5.e12581


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   protocol
   drawing
   data_structures
   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
