
.. include:: links.rst

.. role:: emph

.. _data-structures-doc:

==================================
Data structures and API philosophy
==================================

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

See the `protocol <protocol.html#protocol>`_ for description of how to make these dataframes.

----------------
paralog patterns
----------------

This data structure is how topiary does things like search through descriptions 
of reciprocal BLAST hits to call sequence orthology. This is a dictionary keying
the name of each paralog to either a list of aliases or a compiled regular
expression for the aliases that maps to the name. (If you send in a paralog
pattern as a list of aliases, topiary will automatically compile the regular
expressions). When running the seed-to-alignment pipeline, this structure is
automatically generated from the seed dataframe. 

A quick example. This dictionary indicates that the strings :code:`"MD2"`,
:code:`"ESOP1"`, etc. should all be interpreted as really being :code:`"LY96"`. 

.. code-block:: python
  
  pp = {"LY96":["MD2",
                "myeloid differentiation protein 2",
                "ESOP1",
                "lymphocyte antigen 96"],
        "LY86":["Lymphocyte Antigen 86",
                "Myeloid Differentiation Protein-1",
                "MD1",
                "RP105-associated 3",
                "MMD-1"]
        }

It can be compiled into a set of regular expressions by:

.. code-block:: python

  topiary.io.load_paralog_patterns(pp)
  
  # gives 
  {'LY96': re.compile(r'esop[\ \-_\.]*1|ly[\ \-_\.]*96|lymphocyte[\ \-_\.]*antigen[\ \-_\.]*96|md[\ \-_\.]*2|myeloid[\ \-_\.]*differentiation[\ \-_\.]*protein[\ \-_\.]*2',re.IGNORECASE|re.UNICODE),
   'LY86': re.compile(r'ly[\ \-_\.]*86|lymphocyte[\ \-_\.]*antigen[\ \-_\.]*86|md[\ \-_\.]*1|mmd[\ \-_\.]*1|myeloid[\ \-_\.]*differentiation[\ \-_\.]*protein[\ \-_\.]*1|rp[\ \-_\.]*105[\ \-_\.]*associated[\ \-_\.]*3',re.IGNORECASE|re.UNICODE)}


A couple of notes on how this is compiled. 

1. topiary automatically adds different separators to the regular expression.
   For example, for :code:`"MD2"` topiary will look for :code:`"MD2"`,
   :code:`"MD 2"`, :code:`"MD-2"`, :code:`"MD_2"`, and :code:`"MD.2"`. It inserts
   separators between letters and numbers or any time there is a space/separator in the alias. 
2. topiary does it's best to make sure the regular expressions are unique and 
   will throw an error if the regex for one paralog pulls up a different paralog.
3. topiary puts the name itself (e.g., :code:`"LY96"` above) and adds it into 
   the compiled regular expression. 

If you want to send in your own regular expressions, you can either create 
a paralog pattern dictionary manually (using :code:`re.compile`) or by creating
a compiled paralog pattern dictionary by calling 
`topiary.io.load_paralog_patterns <topiary.io.html#topiary.io.paralog_patterns.load_paralog_patterns>`_
using it's more flexible options. 



API philosophy
==============

---------
Pipelines
---------

The basic idea for the topiary API is to have each function take a topiary
dataframe as an argument and to then return a modified *copy* of that dataframe
as an output. This allows one to write pipelines, as well as easily save out
intermediate steps. The basic flow goes something like:

.. code-block:: python

  df = topiary.do_something(df,args)
  topiary.write_dataframe(df,"current-state.csv")

The main topiary functions take a topiary dataframe as their first argument,
other arguments needed by the function, and then return an appropriately
modified copy of the dataframe. Topiary functions generally modify dataframes by
adding columns with new information and/or by setting the :code:`keep` column to
:code:`True` or :code:`False`. The modified dataframe can then be written out to
a csv file to preserve the current state of the dataframe.

The following code block shows the core of the alignment redundancy reduction
pipeline as one might run it via the API:

.. code-block:: python

  import topiary

  df = topiary.read_dataframe("some_dataframe.csv")

  df = topiary.quality.shrink_in_species(df)
  topiary.write_csv(df,"after-first-shrink.csv")

  df = topiary.quality.shrink_redundant(df)
  topiary.write_csv(df,"after-second-shrink.csv")

  df = topiary.quality.shrink_aligners(df)
  topiary.write_csv(df,"after-third-shrink.csv")


---------------
Run directories
---------------

For the wrapped software that generates output in directories (basically, 
everything after the seed to alignment step), topiary uses a stereotyped 
format for all directories. 

**run_directory**

  + *input*: input files for the calculation
  + *working*: temporary files used when doing the calculation
  + *output*: final output files for the calculation
  + *run_parameters.json*: file holding the run parameters

This is managed by a private API class (:code:`topiary._private.supervisor.Supervisor`)
that can read, write, and manipulate these directories. This class is
documented in the source code and unit tests; advanced users can check those out
if they want to run their own Supervisor-enabled function calls. It remains a private
API class because we do not guarantee it's current functionality in future
release. Any user-written code employing a Supervisor may break without warning. 