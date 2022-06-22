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

Python framework for doing ancestral sequence reconstruction using pandas_
dataframes and ete3_ trees as the primary data structures.

`API reference`_ | :ref:`Installation` | :ref:`Protocol` | :ref:`How to cite` | :ref:`Data Structures`

.. _Installation:

Installation
============

.. note::
  A pip/conda installation protocol is coming soon.


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


.. _Protocol:

Protocol
========

-------------------------
1. Prepare seed dataframe
-------------------------

The first step in a topiary ASR calculation is constructing a
:ref:`seed dataframe`. This table defines the protein family members
(paralogs_) of interest and which species have these proteins (the
taxonomic distribution of the family). topiary will use this to automatically
find and download sequences to put into the alignment. The table can be prepared
in a spreadsheet program (Excel or LibreOffice), a text editor, or
programmatically via pandas_.

An example for the LY86/LY96 protein family is shown below. The full
spreadsheet can be downloaded `here <_static/data/seed-dataframe_example.csv>`_.

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

Protocol for preparing the table:

#. :emph:`Choose the paralogs of interest for your ASR calculation.` In our example,
   we included two paralogs: LY86_ and LY96_, a pair of closely related innate
   immune proteins we study in the lab. The choice of paralogs sets the scope
   of the evolutionary study. In our experience, you'll want ~1-10 paralogs for
   a robust ASR investigation. As you add more paralogs, you need more sequences
   to resolve the evolutionary tree, making the calculation progressively slower.
#. :emph:`Determine the taxonomic distribution of the protein family.` LY86 and LY96 are
   found across bony vertebrates (humans and bony fish, but not sharks). If you
   unsure of the taxonomic distribution of your proteins of interest, you can
   BLAST the protein sequences against the `nr-clustered`_ BLAST database,
   setting the :code:`Max target sequences` parameter to 1000 or more. On the
   resulting output page, the "Taxonomy" tab will show which organisms have
   hits.
#. :emph:`Choose two or three species with well-annotated genomes` that span the whole
   taxonomic distribution of your proteins of interest. For LY86 and LY96, we
   selected humans, chickens, and zebrafish, covering the breadth of species
   over which these proteins are found. (Choosing humans, chimps, and gorillas
   would be a poor choice, as this covers only primates; even choosing humans,
   mice, and chickens would be non-optimal, as this covers only aminotes). The
   NCBI BLAST "Taxonomy" report can be helpful in this regard: select species
   from the highest level of the hierarchy and you'll end up with good
   taxonomic coverage.
#. :emph:`Download sequences for each paralog` from each species and put it into the
   table. Our usual source for these seed sequences is uniprot_. Generally,
   you'll want the canonical_ sequence rather than an isoform. These sequences
   can come from anywhere; they do not even have to be from a database.
#. :emph:`Compile a list of aliases for each paralog`. Annoyingly, the same protein
   can have different names across different databases/species. By using a
   human-curated list of aliases, topiary is more effective at identifying
   sequences that really correspond to the paralogs of interest. When creating
   this list, keep in mind:

   + separate different aliases with `;`
   + these aliases are case-insensitive (:code:`MD2`, :code:`md2`, and
     :code:`Md2` are the same)
   + topiary automatically tries different separators. For example, for :code:`MD2`
     topiary will look for :code:`MD2`, :code:`MD 2`, :code:`MD-2`, :code:`MD_2`,
     and :code:`MD.2`. It inserts separators between letters and numbers or any
     time there is a space/separator in the pattern in the alias.
   + make sure to include both the abbreviated and full-length versions of each
     alias (i.e. :code:`MD2` and :code:`myeloid differentiation protein 2`).

   To find aliases, check out the `Also known as`_ field for the gene of interest
   on NCBI, the `Protein names`_ section of the protein's uniprot entry, and/or
   (for proteins found in humans) the relevant genecards_ entry.
#. You can put other information about the sequences (accession, citations, etc.)
   as their own columns in the table. topiary will ignore, but keep, those
   columns.

------------------------
2. Generate an alignment
------------------------

 .. note::
   In the final version of this protocol, we will provide a full example
   jupyter notebook *and* describe how to do the following steps using a
   command line.

In a jupyter notebook, run the
`seed_to_alignment <api/topiary/pipeline/seed_to_alignment.html>`_ function.

.. code-block:: python

  import topiary
  df = topiary.seed_to_alignment("seed-dataframe_example.csv",
                                 out_dir="seed_to_alignment")

.. note::
  This function takes ~5 minutes to run on an M1 macbook pro using the example
  :code:`seed-dataframe_example.csv` file. The timing for the initial NCBI BLAST
  step depends on server capacity and may take awhile. You might want to grab
  a cup of coffee while you wait.

The `seed_to_alignment <api/topiary/pipeline/seed_to_alignment.html>`_ function
will:

+ :emph:`Find paralogs from other species` by using the seed sequences to as BLAST queries
  against the `nr <https://www.ncbi.nlm.nih.gov/refseq/about/nonredundantproteins/>`_
  database.
+ :emph:`Perform initial orthology calls` using
  `reciprocal BLAST <https://www.flyrnai.org/RNAi_orthology.html>`_ against the
  proteomes of the species from the seed dataframe.
+ :emph:`Remove difficult-to-align` sequences from the dataset by performing a
  preliminary alignment using muscle's `super5 <https://drive5.com/muscle5/manual/cmd_super5.html>`_ algorithm,
  then checking alignment quality for each sequence
+ :emph:`Create a reasonably sized alignment` by culling sequences in a
  taxonomically-informed way. Instead of removing sequences based on their
  relative sequence identity, topiary removes taxonomically redundant sequences.
  (For example, if choosing which to remove from a sequences taken from human,
  chimp, or turkey, it would remove the chimp or human because they are more
  closely related to each other than to the turkey).
+ :emph:`Perform a final alignment` using muscle. This alignment is written out
  to :code:`05_alignment.fasta`.

-----------------------------------------------------
3. Visually inspect and (possibly) edit the alignment
-----------------------------------------------------

Check out the alignment in `aliview <https://ormbunkar.se/aliview/>`_.

If you edited the alignment, you need to load the new alignment back into the
topiary dataframe. In a jupyter notebook, run the following:

.. code-block:: python

  import topiary
  df = topiary.load_dataframe("seed_to_alignment/04_aligned-dataframe.csv")
  df = topiary.read_fasta_into(df,EDITED_FASTA_FILE)
  topiary.write_dataframe(df,"07_final-alignment.csv")

------------------------------------------------
4. Infer phylogenetic model, tree, and ancestors
------------------------------------------------

.. note::
  We highly recommend running the following steps on a computing cluster. To
  prepare the computing environment, please follow the installation steps above on
  the cluster.

.. note::
  The following code block will be a single command line call when this is done.

This program will:

+ :emph:`Choose a phylogenetic model` that maximizes the likelihood of observing
  the sequences in your alignment. (Uses RAxML-NG).
+ :emph:`Generate a maximum-likelihood phylogenetic tree` with bootstrap
  supports. (Uses RAxML-NG).
+ :emph:`Reconcile the gene and species trees`. (Uses GeneRax).
+ :emph:`Estimate ancestral sequences` using an empirical Bayes method. (Uses
  RAxML-NG).

The final output will be an `ancestors` directory.

Run the following python script on the cluster using your job scheduler
(slurm, etc.). Change NUM_THREADS_YOU_CHOOSE to the same number of threads you
allocate to the job. This takes ?? hours on our cluster.

.. code-block:: python

  import topiary

  num_threads = NUM_THREADS_YOU_CHOOSE
  df_file = "07_final-alignment.csv"

  topiary.find_best_model(final_df,
                          output="find-model",
                          threads=num_threads)
  topiary.generate_ml_tree(previous_dir="find-model",
                           output="ml-tree",
                           threads=num_threads,
                           bootstrap=True)
  topiary.reconcile(previous_dir="ml-tree",
                    output="reconciled",
                    overwrite=True)
  topiary.generate_ancestors(previous_dir="reconciled",
                             output="ancestors",
                             threads=num_threads,
                             overwrite=True)

--------------------------------------------------------
5. Analyze results and run bootstraps on reconciled tree
--------------------------------------------------------

TODO.




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


.. _Data Structures:

Data structures
===============

Under the hood, topiary uses pandas_ dataframes to manage the phylogenetic data
in the project. For those unfamiliar with dataframes, these are essentially
spreadsheets with a row for each sequence and columns holding various features
of that sequence. These dataframes can be readily written out and read from
spreadsheet files (.csv, .tsv, .xlsx).

There are two key dataframes you will encounter using topiary:

+ *topiary dataframe*: The main structure for holding sequences and information
  about those sequences for the project. You can think of it as a giant
  spreadsheet of sequences.
+ *seed dataframe*: A manually constructed dataframe containing seed sequences
  that topiary will then use to construct a full topiary dataframe for the
  project.

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

Other user-specified columns are allowed. In addition, specific topiary analyses
may add new columns. For example, `recip_blast` will add multiple columns
including `recip_paralog` and `recip_prob_match`.

Reading and writing
-------------------

Topiary dataframes are standard pandas_ dataframes and can thus be written to
and read from various spreadsheet formats. We recommend using topiary's
built-in functions to read and write the dataframes (`topiary.read_dataframe`
and `topiary.write_dataframe`). These functions will preserve/check column
formats etc.

You can manually edit a topiary dataframe using pandas_ operations or using a
spreadsheet program (i.e. Excel). If you manually edit a dataframe, make sure
that all sequences have unique `uid` and that all sequences in the `alignment`
column, if present, have identical length.

Constructing
------------

There are a number of ways to construct a topiary dataframe:

+ :code:`io.df_from_seed`: construct topiary dataframe from a seed dataframe.
+ :code:`io.df_from_blast_xml`: construct topiary dataframe from xml files
  generated by NCBI blast.


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

.. _opentreeoflife: https://tree.opentreeoflife.org/
.. _ott770315: https://tree.opentreeoflife.org/opentree/argus/ottol@770315/Homo-sapiens
.. _ott276534: https://tree.opentreeoflife.org/opentree/argus/ottol@276534/Thermus-thermophilus
.. _pandas: https://pandas.pydata.org/docs/
.. _ete3: http://etetoolkit.org
.. _GeneRax: https://github.com/BenoitMorel/GeneRax
.. _`Open Tree of Life`: https://tree.opentreeoflife.org/
.. _PastML: https://pastml.pasteur.fr
.. _`RAxML-NG`: https://github.com/amkozlov/raxml-ng
.. _Muscle: https://www.drive5.com/muscle/
.. _paralogs: https://en.wikipedia.org/wiki/Sequence_homology#Paralogy
.. _LY96: https://www.uniprot.org/uniprot/Q9Y6Y9
.. _LY86: https://www.uniprot.org/uniprot/O95711
.. _`NCBI BLAST`: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins
.. _`nr-clustered`: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
.. _uniprot: https://www.uniprot.org
.. _canonical: https://www.uniprot.org/help/canonical_and_isoforms
.. _genecards: https://www.genecards.org
.. _`Also known as`: https://www.ncbi.nlm.nih.gov/gene/23643
.. _`Protein names`: https://www.uniprot.org/help/protein_names
.. _`API reference` : topiary.html


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   topiary


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
