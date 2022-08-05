
.. include:: links.rst

============
Installation
============

Requirements
============

+ Core scientific python libraries:

  + `Python >= 3.8 <https://www.python.org/downloads/>`_
  + `numpy <https://numpy.org>`_
  + pandas_
  + matplotlib_

+ Tree manipulation/drawing:

  + `ete3 <http://etetoolkit.org/download/>`_
  + `toytree <https://toytree.readthedocs.io/en/latest/3-installation.html>`_

+ Packages used for tree/ancestor inferences:

  + `NCBI BLAST+ >= 2.2 <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_
  + `muscle >= 5.0 <https://github.com/rcedgar/muscle/releases/>`_
  + `GeneRax >= 2.0 <https://github.com/BenoitMorel/GeneRax/releases/>`_
  + `RAxML-NG >= 1.1 <https://github.com/amkozlov/raxml-ng/releases/>`_
  + `pastml >= 1.9 <https://pastml.pasteur.fr>`_


macOS and linux
===============

We recommend installing topiary via conda_. If you do not have conda installed,
download and install miniconda_ before proceeding.

To prevent interference with other packages, we recommend installing topiary in
its own conda environment. To do so, open a terminal and type:

.. code-block:: shell-session

  conda create -n topiary topiary-asr

This will automatically install all of the libraries and software to run topiary.

You can make sure that all of the pipeline dependencies installed correctly by
opening a terminal and typing:

.. code-block:: shell-session

  conda activate topiary
  topiary-test-installed

This will spit out a list of software that topiary detected.

If you want to use topiary interactively in jupyter notebooks, you will need to
install jupyter-lab. In a terminal, type:

.. code-block:: shell-session

  conda activate topiary
  conda install jupyter-lab

-------
Windows
-------

We recommend installing topiary via conda_. If you do not have conda installed,
download and install miniconda_ before proceeding.

To prevent interference with other packages, we recommend installing topiary in
its own conda environment. To do so, open the conda terminal and type:

.. code-block:: shell-session

  conda create -n topiary topiary-asr
