
.. include:: links.rst

.. role:: emph

.. _data-structures-doc:

============================
Patterns and data structures
============================

Patterns
===============

Operates on a copy


.. code::

  df = topiary.do_something(df,args)





Data structures
===============

Under the hood, topiary uses `pandas <pandas-link_>`_ dataframes to manage the
phylogenetic data in the project. For those unfamiliar with dataframes, these are
essentially spreadsheets with a row for each sequence and columns holding various
features of that sequence. These dataframes can be readily written out and read
from spreadsheet files (.csv, .tsv, .xlsx).

Topiary is built around two types of dataframes:

+ :emph:`seed dataframe`: A manually constructed dataframe containing seed
  sequences that topiary uses as input to construct a full topiary dataframe for
  the project.
+ :emph:`topiary dataframe`: The main structure for holding sequences and
  information about those sequences. Each step in the pipeline edits, saves
  out, and then returns the main dataframe. This allows one to follow
  the steps and/or manually introduce changes.

-----------------
topiary dataframe
-----------------

A topiary dataframe must have three columns:

+ :code:`name`: a name for the sequence. This does not have to be unique.
+ :code:`sequence`: the amino acid sequence. This does not have to be unique.
+ :code:`species`: the species name for this sequence (binomial, i.e. *Homo sapiens*
  or *Thermus thermophilus*).

Topiary will automatically add a few more columns if not present.

+ :code:`keep`: a boolean (True/False) column indicating whether or not to
  use the sequence in the analysis. Topiary will not delete a sequence from the
  dataset, but instead set :code:`keep = False`.
+ :code:`uid`: a unique 10-letter identifier for this sequence.

.. danger:: uid values should never be modified by the user.

+ :code:`ott`: The opentreeoflife_ reference taxonomy identifier for the
  sequence species. This will have the form ottINTEGER (i.e. ott770315_ for
  *Homo sapiens* and ott276534_ for *Thermus thermophilus*).

Topiary reserves a few more columns that may or may not be used:

+ :code:`alignment`: an aligned version of the sequence. All sequences in the
  alignment column must have the same length.
+ :code:`always_keep`: a boolean (True/False) column indicating whether or not
  topiary can drop the sequence from the analysis.

In addition, specific topiary analyses may add new columns. For example,
:code:`recip_blast` will add multiple columns such as :code:`recip_paralog` and
:code:`recip_prob_match`.

Other user-specified columns are allowed.


Constructing
------------

There are two basic ways to construct a topiary dataframe:

+ :code:`io.df_from_seed`: construct topiary dataframe from a seed dataframe.
  Depending on the options selected, topiary will add sequences using BLAST or
  will read sequences from a list of pre-prepared BLAST xml files.
+ Construct the dataframe manually.


Reading and writing
-------------------

Topiary dataframes are standard `pandas <pandas-link_>`_ dataframes and can thus be written to
and read from various spreadsheet formats. We recommend using topiary's
built-in functions to read and write the dataframes (`topiary.read_dataframe`
and `topiary.write_dataframe`). These functions will preserve/check column
formats etc.


Editing
-------

You can manually edit a topiary dataframe using `pandas <pandas-link_>`_ operations or using a
spreadsheet program (i.e. Excel). If you manually edit a dataframe, make sure
that all sequences have unique `uid` and that all sequences in the `alignment`
column, if present, have identical length.


.. _seed dataframe:

--------------
seed dataframe
--------------

A seed dataframe must have four columns:

+ :code:`name`: name of each sequence. This will usually be a short, useful
  name for the paralog.
+ :code:`species`: species names for seed sequences in binomial format (i.e.
  *Homo sapiens* or *Thermus thermophilus*).
+ :code:`aliases`: other names for each protein that may be used in various
  databases/species, separated by :code:`;`.
+ :code:`sequence`: amino acid sequences for these proteins.


Example seed dataframe
----------------------

+------+-------------------------------------------------------------------+--------------+------------+
| name | aliases                                                           | species      | sequence   |
+------+-------------------------------------------------------------------+--------------+------------+
| LY96 | lymphocyte antigen 96;MD2;ESOP1;Myeloid Differentiation Protein-2 | Homo sapiens | MLPFLFF... |
+------+-------------------------------------------------------------------+--------------+------------+
| LY96 | lymphocyte antigen 96;MD2;ESOP1;Myeloid Differentiation Protein-2 | Danio rerio  | MALWCPS... |
+------+-------------------------------------------------------------------+--------------+------------+
| LY86 | lymphocyte antigen 86;MD1;Myeloid Differentiation Protein-1       | Homo sapiens | MKGFTAT... |
+------+-------------------------------------------------------------------+--------------+------------+
| LY86 | lymphocyte antigen 86;MD1;Myeloid Differentiation Protein-1       | Danio rerio  | MKTYFNM... |
+------+-------------------------------------------------------------------+--------------+------------+

See protocol for description of how to make these dataframes.


Run directories
===============

run_directory
+ input
+ working
+ output
+ run_parameters.json

.. code_block::
    run_parameters = {"files":["files"] <<- relative
                      "version":
                      "calc_type":
                      "kwargs":
                      "start_time":
                      "end_time":
                      "model":
                      "cmd":
                      "previous_entries":[]}
