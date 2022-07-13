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

------------
Requirements
------------

+ `Python >= 3.8 <https://www.python.org/downloads/>`_
+ pandas_
+ ete3_
+ `numpy <https://numpy.org>`_
+ `NCBI BLAST+ >= 2.2 <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/>`_
+ `GeneRax >= 2.0 <https://github.com/BenoitMorel/GeneRax/releases/>`_
+ `RAxML-NG >= 1.1 <https://github.com/amkozlov/raxml-ng/releases/>`_
+ `muscle >= 5.0 <https://github.com/rcedgar/muscle/releases/>`_

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


.. _Protocol:

Protocol
========

There are some steps in phylogenetics that are best done by humans; most steps
are best done by computers. Our philosophy is to facilitate the human-centric
steps, then automate the rest.

The most important task in phylogenetics is defining the problem. What
ancestors do you want to reconstruct? What modern proteins are best
characterized and most relevant to interpreting results with ancestors?



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

Protocol for preparing the table
--------------------------------

#. :emph:`Choose the paralogs of interest for your ASR calculation.` In our example,
   we included two paralogs: LY86_ and LY96_, a pair of closely related innate
   immune proteins we study in the lab. The choice of paralogs sets the scope
   of the evolutionary study. In our experience, you'll want ~1-10 paralog
   sequences in your seed dataframe for a robust ASR investigation. As you add
   more paralogs, you need more sequences to resolve the evolutionary tree,
   making the calculation progressively slower.
#. :emph:`Determine the taxonomic distribution of the protein family.` LY86 and LY96 are
   found across bony vertebrates (humans and bony fish, but not sharks). If you are
   unsure of the taxonomic distribution of your proteins of interest:

   + BLAST known protein sequences against the `nr-clustered`_ BLAST database,
     setting the :code:`Max target sequences` parameter to 1000 or more. On the
     resulting output page, the "Taxonomy" tab will show which organisms have
     hits.
   + Take a few representative sequences from the most divergent species in the
     outputs. BLAST those against the NCBI `nr`_ database...

#. :emph:`Choose two or three species with well-annotated genomes` that span the
   whole taxonomic distribution of your proteins of interest. For LY86 and LY96,
   we selected humans, chickens, and zebrafish, covering the breadth of species
   over which these proteins are found. (Choosing humans, chimps, and gorillas
   would be a poor choice, as this covers only primates; even choosing humans,
   mice, and chickens would be non-optimal, as this covers only amniotes). The
   NCBI BLAST "Taxonomy" report can be helpful in this regard: select species
   from the highest level of the hierarchy and you'll end up with good
   taxonomic coverage.
#. :emph:`Download sequences for each paralog` from each species and put it into
   the table. Our usual source for these seed sequences is uniprot_. Generally,
   you'll want the canonical_ sequence rather than an isoform. These sequences
   can come from anywhere; they do not even have to be from a database.
#. :emph:`Compile a list of aliases for each paralog`. Annoyingly, the same protein
   can have different names across different databases/species. By using a
   human-curated list of aliases, topiary is more effective at identifying
   sequences that really correspond to the paralogs of interest. When creating
   this list, keep in mind:

   + separate different aliases with `;`
   + these aliases are case-insensitive (i.e. :code:`MD2`, :code:`md2`, and
     :code:`Md2` are equivalent)
   + topiary automatically tries different separators. For example, for :code:`MD2`
     topiary will look for :code:`MD2`, :code:`MD 2`, :code:`MD-2`, :code:`MD_2`,
     and :code:`MD.2`. It inserts separators between letters and numbers or any
     time there is a space/separator in the pattern in the alias.
   + make sure to include both the abbreviated and full-length versions of each
     alias (i.e. :code:`MD2` and :code:`myeloid differentiation protein 2`).

   To find aliases, you can check out the `Also known as`_ field for the gene of
   interest on NCBI, the `Protein names`_ section of the protein's uniprot
   entry, a genecards_ entry (for proteins found in humans), and/or primary
   literature.
#. You can put other information about the sequences (accession, citations, etc.)
   as their own columns in the table. topiary will ignore, but keep, those
   columns.

--------------------
2. Seed to alignment
--------------------

Generate an alignment on the command line.

Code
----

.. code-block:: shell-session

  topiary-seed-to-alignment seed-dataframe_example.csv output

.. note::

    The timing for the initial NCBI BLAST step depends on server capacity and
    may take awhile. If the search hangs for a long time (say, 20 min), you can
    stop the program (hit ⌘+C or CTRL+C) and try again. Unfortunately, the
    time required for this search is outside topiary's control. For ~10 seeds
    generating an alignment of 500 sequences, the topiary code itself takes ~5
    minutes to run on a 2021 M1 Macbook Pro.


.. note::

  If this function returns many fewer sequences than your target alignment
  size, the output will provide information that can help you understand why.
  One common problem few reciprocal blast hits because you have not specified
  all of the aliases for your protein in the NCBI database. To fix this
  problem, open the file :code:`output/02_recip-blast-dataframe.csv` and look in
  the *recip_found* and *recip_hit* columns. If you find many entries that
  have :code:``recip_found == False`` but look like they should match in the
  *recip_hit* column, you can add a new alias that would match these missed hits
  to the *aliases* column in your seed dataframe.

Details
-------

This function does the following:

+ :emph:`Find paralogs from other species` using the seed sequences as BLAST queries
  against the `nr <https://www.ncbi.nlm.nih.gov/refseq/about/nonredundantproteins/>`_
  database.
+ :emph:`Perform initial orthology calls` using
  `reciprocal BLAST <https://www.flyrnai.org/RNAi_orthology.html>`_ against the
  proteomes of the species from the seed dataframe.
+ :emph:`Download a species tree for sequence hits` from the opentreeoflife_
  taxonomic database.
+ :emph:`Create rough alignments of orthologous sequences` muscle's
  `super5 <https://drive5.com/muscle5/manual/cmd_super5.html>`_ algorithm.
  Topiary uses this information when lowering sequence redundancy on the next
  step.
+ :emph:`Create a reasonably sized alignment` by culling sequences in a
  taxonomically-informed way. Instead of removing sequences based on their
  relative sequence identity, topiary removes taxonomically redundant sequences
  using the species tree as a guide. For example, if choosing one sequence to
  remove from a set of sequences taken from human, chimp, or turkey, it would
  remove the chimp or human because they are more closely related to each other
  than to the turkey. It would select between the human and chimp sequences
  based their relative alignment quality in alignment from the last step.
+ :emph:`Perform a final alignment` using muscle. This alignment is written out
  to :code:`05_alignment.fasta`.

Options
-------

To find out and set the options for this pipeline, use the :code:`--help` flag.
Running the following:

.. code-block:: shell-session

  topiary-seed-to-alignment --help


Will return information about the flags that can be specified in the program.
You can access the same information in the
`topiary documentation <topiary.pipeline.html#module-topiary.pipeline.seed_to_alignment>`_

Output
------

This will output a directory with the following files:

+ 00_xx. A copy of the seed sequence file.
+ 01_initial-dataframe.csv. All sequences downloaded from NCBI.
+ 02_recip-blast-dataframe.csv. Results of reciprocal BLAST.
+ 03_sampled-dataframe.csv. Dataframe with sequence redundancy lowered.
+ 04_aligned-dataframe.csv. Dataframe with final alignment.
+ 05_alignment.fasta. Dataframe alignment written out to a fasta file.

There are other files in this directory (.faa.gz and blast_db.* files). These
were used for the reciprocal BLAST calculation and may be deleted if desired.


-----------------------------------------------------
3. Visually inspect and (possibly) edit the alignment
-----------------------------------------------------

Alignments aren't always perfect. So we can edit the alignment.

Check out the alignment in `aliview <https://ormbunkar.se/aliview/>`_.

Once you have edited the alignment, you need to load the new alignment back into the
topiary dataframe.

Command line instructions
-------------------------

.. code-block:: shell-session

  topiary-read-fasta-into 04_aligned-dataframe.csv EDITED_FASTA NEW_CSV


Jupyter instructions
--------------------

In a jupyter notebook, run the following:

.. code-block:: python

  import topiary
  df = topiary.load_dataframe(LAST_CSV_FILE)
  df = topiary.read_fasta_into(df,EDITED_FASTA_FILE)
  topiary.write_dataframe(df,NEW_CSV_FILE)

------------------------------------------------
4. Infer phylogenetic model, tree, and ancestors
------------------------------------------------

.. note::
  We highly recommend running the following steps on a computing cluster. To
  prepare the computing environment, please follow the installation steps above on
  the cluster.

Assuming you are running this on a computing cluster, you'll need to specify
the resources available for the calculation.

.. code-block:: shell-session

  #!/bin/bash -l
  #SBATCH --account=harmslab
  #SBATCH --job-name=topiary
  #SBATCH --output=hostname.out
  #SBATCH --error=hostname.err
  #SBATCH --partition=long
  #SBATCH --time=07-00:00:00
  #SBATCH --nodes=1
  #SBATCH --ntasks-per-node=1
  #SBATCH --cpus-per-task=28

  ~/topiary/bin/topiary-alignment-to-ancestors 04_aligned-dataframe.csv --out_dir output --num_threads 28

.. code-block:: shell-session

  qsub launch_topiary.srun


Details
-------

This pipeline will:

+ :emph:`Choose a phylogenetic model` that maximizes the likelihood of observing
  the sequences in your alignment. Topiary generates a rough tree tree
  (maximum parsimony) and then optimizes the tree branch lengths using ~300
  different phylogenetic models implemented in raxml-ng. It then selects between
  those models using an AIC test.
+ :emph:`Generate a maximum-likelihood phylogenetic tree` with bootstrap
  supports. (Uses RAxML-NG).
+ :emph:`Reconcile the gene and species trees`. (Uses GeneRax).
+ :emph:`Estimate ancestral sequences` using an empirical Bayes method. (Uses
  RAxML-NG).

The final output will be an `03_ancestors/output` directory.

Output
------

00_find-model  01_ml-tree  02_reconciliation  03_ancestors

Each of these will have an *output* directory.

+ dataframe.csv
+ run_parameters.json
+ 00_find_model

  - model-comparison.csv. Has all model comparisons, with the selected model at
    the top of the file. The right-most column (:code:`p`) is the probability
    that this model minimizes the information loss relative to the other tested
    models. A value of 1.0 implies high confidence in the selected model.

+ 01_ml. Maximum likelihood tree estimation.

  - tree.newick. Maximum likelihood tree with branch lengths in machine-readable
    format. If bootstrap replicates were done, this tree will also have branch
    supports.
  - summary-tree.pdf. Maximum likelihood tree as a vector-graphics image. If
    bootstrap replicates were performed, this tree will have bootstrap support
    encoded by color. Note that this tree is not rooted.

+ 02_reconciliation. Gene/species tree reconcilation results.

  - tree.newick. Reconciled gene/species tree with branch lengths in
    machine-readable format. If bootstrap replicates were done, this tree will
    also have branch supports.
  - summary-tree.pdf. Reconciled gene/species tree as a vector-graphics image.
    Evolutionary events (duplication, speciation, etc.) are labeled on the tree.
    Note that this tree *is* rooted. If bootstrap replicates were performed,
    this tree will have bootstrap support encoded by color.
  - reconcilations. Directory with output from the reconcilations calculation.

+ 03_ancestors. Reconstructed ancestral proteins.


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
such as `recip_paralog` and `recip_prob_match`.

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
.. _`miniconda` : https://docs.conda.io/en/latest/miniconda.html


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   topiary


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
