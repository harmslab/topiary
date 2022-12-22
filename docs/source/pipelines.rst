
.. role:: emph

.. role:: raw-html(raw)
    :format: html

.. include:: links.rst

.. _pipelines-doc:

================
Pipeline details
================

This page provides more details about the steps and intermediate outputs of the
pipelines described in the topiary :ref:`protocol<protocol-doc>`. 

.. _align-script-details:

Seed to alignment
=================

Seed to alignment starts with a set of seed sequences and returns a high quality
initial alignment. The steps it does are:

#. Use the seed sequences as BLAST queries against relevant databases. (Default
   is the NCBI *nr* database).
#. Discard sequences that do not return a seed sequence when used as reciprocal
   BLAST queries against proteomes from key species.
#. For non-microbial datasets, remove sequences from species that cannot be
   resolved on the most recent Open Tree of Life synthetic tree.
#. Remove similar sequences within each species using a strict sequence identity
   cutoff (default = 0.95). This removes isoforms and recent lineage-specific
   duplications.
#. Calculate a target alignment size based on the length of the longest seed 
   sequence. By default, aim to build a dataset with one sequence per amino acid
   in that sequence. Then multiply this by 1.1 so we can remove sequences after
   budgeting and still have an alignment of the desired size.
#. For a microbial dataset, identify a sequence identity cutoff that yields an
   alignment of the right size. For a non-microbial dataset, identify blocks of
   sequences to merge that sample the species tree. (See the `topiary paper <topiary-link_>`_
   for details). 
#. Align all sequences in the current dataset, which will have ~1.1 times the
   target alignment size. Remove the worst aligning sequences. This is done
   + Remove the sequences with the most characters in non-dense columns (drop worst 2.5%).
   + Remove the sequences with the most missing dense columns (drop worst 2.5%).
#. Re-align the dataset. It should now have ~1.05 times the target alignment
   size and be ready for manual alignment editing. 


The pipeline does writes out five *.csv* files over the course of this analysis,
allowing one to track what changes are made. Topiary will add sequences and/or
columns at each step. Until the final step, it does not delete sequences, but
rather sets the :code:`keep` column to :code:`False` when a sequence is removed.

+ :emph:`Finds paralogs from other species` using the seed sequences as BLAST queries
  against the `NCBI non-redundant database <blast-nr_>`_. The taxonomic scope is
  defined by the key species in the seed dataset. For bacterial datasets, the 
  taxonomic scope will be all bacteria; for archaeal datasets, the taxonomic scope
  will be all archaea. 

  + Output: *01_initial-dataframe.csv*
  + Function: `topiary.ncbi.ncbi_blast <topiary.ncbi.blast.html#topiary.ncbi.blast.ncbi.ncbi_blast>`_

+ :emph:`Makes orthology calls` using `reciprocal BLAST <https://www.flyrnai.org/RNAi_orthology.html>`_
  against the proteomes of key species from the seed dataset.

  + Output: *02_recip-blast-dataframe.csv*
  + Function: `topiary.ncbi.recip_blast <topiary.ncbi.blast.html#topiary.ncbi.blast.recip.recip_blast>`_

+ :emph:`Decreases the size of the dataset` taking into account taxonomic sampling,
  quality, and alignment score of all sequences. For non-microbial datasets, 
  topiary uses the species tree to select phylogenetically informative sequences;
  for microbial datasets, topiary lowers redundany based on sequence similarity
  alone. 

  + Output: *03_shrunk-dataframe.csv*
  + Function: `topiary.shrink_dataset <topiary.quality.html#topiary.quality.shrink.shrink_dataset>`_

+ :emph:`Generates a draft alignment` using `Muscle5 <muscle-link_>`_.

  + Output: *04_aligned-dataframe.csv*
  + Function: `topiary.align <topiary.muscle.html#topiary.muscle.muscle.align>`_

+ :emph:`Polishes alignment` by removing worst aligning sequences.

  + Output: *05_clean-aligned-dataframe.csv*
  + Function: `topiary.quality.polish_alignment <topiary.quality.html#topiary.quality.polish.polish_alignment>`_

+ :emph:`Writes alignment in fasta format` for loading into an alignment viewer.

  + Output: *06_alignment.fasta*
  + Function: `topiary.write_fasta <topiary.io.html#topiary.io.alignments.write_fasta>`_



.. _ali-to-anc-pipeline:

Alignment to Ancestors
======================

This pipeline does the following steps. 

+ Infer the maximum likelihood model of sequence evolution
+ Infer a maximum likelihood gene tree
+ Infer ancestors on the maximum likelihood gene tree
+ Reconcile the gene and species trees (for non-microbial proteins)
+ Generate bootstrap replicates for the maximum likelihood gene tree

For each step, we describe what the software does, as well as showing example
output files.

.. note::

  The *.newick* tree files are copied into the :code:`output` directory from
  all previous steps. (This means, for example, that the :code:`output`
  directory from the last step will have all *.newick* files generated
  previously.) In what follows, we only describe the *new* output that
  will be present after each step.


+ :emph:`Choose a phylogenetic model` that maximizes the likelihood of observing
  the sequences in your alignment. Topiary uses RAxML-NG to generate a maximum 
  parsimony tree and then optimizes the tree branch lengths using 360
  different phylogenetic models implemented in RAxML-NG. It selects between
  those models using a corrected AIC test. The following table is an example of
  the resulting output.

  +---+-----------+---------+---------+----------+---------+---+---------+
  |   |model      |L        |AIC      |AICc      |BIC      |N  |p        |
  +---+-----------+---------+---------+----------+---------+---+---------+
  |189|LG+G8      |-86776.11|175548.23|2169552.23|179924.65|998|1.00     |
  +---+-----------+---------+---------+----------+---------+---+---------+
  |351|WAG+G8     |-86998.12|175992.25|2169996.25|180368.67|998|3.82e-97 |
  +---+-----------+---------+---------+----------+---------+---+---------+
  |333|VT+G8      |-87026.31|176048.63|2170052.63|180425.06|998|2.17e-109|
  +---+-----------+---------+---------+----------+---------+---+---------+
  |360|LG4M       |-87163.38|176322.77|2170326.77|180699.20|998|6.44e-169|
  +---+-----------+---------+---------+----------+---------+---+---------+
  |63 |DEN+G8     |-87215.32|176426.64|2170430.64|180803.07|998|1.79e-191|
  +---+-----------+---------+---------+----------+---------+---+---------+
  |81 |Blosum62+G8|-87273.60|176543.21|2170547.21|180919.64|998|8.74e-217|
  +---+-----------+---------+---------+----------+---------+---+---------+
  |...|...        |...      |...      |...       |...      |...|...      |
  +---+-----------+---------+---------+----------+---------+---+---------+

  + Output: :code:`00_find-model/output/`

    + *model_comparison.csv*: This file has models ranked form best to worst.
      The best model has the highest value in the probability (:code:`p`) column.
      For each model, topiary calculates :math:`w_{i} = exp(({AICc}_{i} - {AICc}_{best})/2)`.
      The probability is :math:`P_{i} = w_{i} / \sum _{j=0}^{j < M} w_{j}`, where
      :math:`M` is the number of models. Other statistics are reported as well,
      including the likelihood (:code:`L`), raw AIC, BIC, and number of model
      parameters (:code:`N`).
    + *dataframe.csv*: Current topiary dataframe.

  + Function: `topiary.find_best_model <topiary.raxml.html#topiary.raxml.model.find_best_model>`_

+ :emph:`Generate a maximum-likelihood gene tree`. Uses RAxML-NG to find the
  maximum likelihood gene tree using the model determined in the previous
  step. It starts the inference from ten different parsimony trees and ten random
  trees, then optimizes the trees using default RAxML-NG moves (NR-FAST and SPR).
  This tree an example of "summary-tree.pdf" generated by this step.

  .. image:: _static/img/ali-to-anc/01_gene-tree.svg
    :align: center
    :alt: maximum likelihood gene tree
    :width: 75%

  + Output: :code:`01_gene-tree/output/`

    + *summary-tree.pdf*: Drawing of ML gene tree with branch lengths. Tree is
      drawn with midpoint rooting.
    + *gene-tree.newick*: ML gene tree with branch lengths in newick format.
      Tips are topiary UIDs. Internal nodes are unlabeled.
    + *dataframe.csv*: Current topiary dataframe.

  + Function: `topiary.generate_ml_tree <topiary.raxml.html#topiary.raxml.tree.generate_ml_tree>`_

+ :emph:`Infers ancestral sequences` on the gene tree using RAxML-NG and the
  maximum likelihood substitution model using RAxML-NG. Gaps are inferred using
  parsimony via PastML. This tree an example of "summary-tree.pdf" generated by
  this step. The tree is rooted by the midpoint method. The graph shows an
  example ancestor sequence summary. The text shows the corresponding
  reconstructed ancestor and its altAll versions as a fasta file.

  .. image:: _static/img/ali-to-anc/02_gene-tree-ancestors.svg
    :align: center
    :alt: ancestors drawn on gene tree
    :width: 75%

  :raw-html:`<br />`

  .. image:: _static/img/ali-to-anc/02_gene-tree-ancestors_plot.svg
    :align: center
    :alt: Graph showing ancestral quality
    :width: 85%

  :raw-html:`<br />`

  .. code-block::

    >anc9|avgPP:0.782|lnPP:-65.837|num_ambig:30|num_ambig_gaps:0
    MKVFFTLLFVLILF-CS----G------ESKEWPTHTICNTSDLEVYYKSCDPL--Q-DVGVSISPCSKSMTEDINVRVALLLRQDIKELYLNLDLYIN----GLHVLSY--D--YPLCEPS-F-PRFTFCGRRKGELISFEGPVKSNIQTIPKGEY-NVSLELFNEDN--Y--TIACANITLI---------------S---R
    >anc9_altAll|avgPP:0.760|lnPP:-76.920|num_ambig:30|num_ambig_gaps:0
    MKIYFTLLLVLILF-CT----G------ESKEWPTHTLCNTSNLEVYYRSCDPL--Q-DIGLSISPCSKSLTEDINIRIALILRQNINELYLNIDVFIN----GLKVLNY--D--YPLCEPS-F-PKFTFCGRKKGEMISFEGPIKSNVQTLPKGEF-NVTLELFNEDN--Y--TIACANVTLI---------------N---R

  + Output: :code:`02_gene-tree-ancestors`

    + *summary-tree.pdf*: Drawing of gene tree with branch lengths and ancestors
      with posterior probabilities as a color map.
    + *gene-tree_anc-pp.newick*: Rooted, reconciled tree with branch
      lengths in newick format. Tips are topiary UIDs. Internal nodes are
      ancestral posterior probabilities.
    + *gene-tree_anc-label.newick*: Rooted, reconciled tree with branch
      lengths in newick format. Tips are topiary UIDs. Internal nodes are
      ancestor names.
    + *gene-tree_ancestors/*: The *ancX.pdf* files are pdf summaries of
      ancestral reconstructions for each ancestor as labeled in *summary-tree.pdf*.
      An example is shown above. These graphs show the support for: the most
      likely reconstruction at each site (black circles), the next most likely
      reconstruction (red circles), the locations of gaps (gray shading), and
      the location of ambiguous gaps (purple dashes). Statistics at the top
      indicate ancestor quality.

  + Function: `topiary.generate_ancestors <topiary.raxml.html#topiary.raxml.ancestors.generate_ancestors>`_

+ :emph:`Reconciles the gene and species trees`. Uses GeneRax to reconcile the
  gene and species trees. Uses default GeneRax SPR moves. User can set
  whether or not to allow horizontal gene transfer (the :code:`UndatedDTL` or
  :code:`UndatedDL` GeneRax models, respectively). This tree an example of
  "summary-tree.pdf" generated by this step.

  .. image:: _static/img/ali-to-anc/03_reconciled-tree.svg
    :align: center
    :alt: maximum likelihood reconciled gene/species tree
    :width: 75%

  + Output: :code:`03_reconciled-tree/output/`

    + *summary-tree.pdf*: Drawing of reconciled tree with branch lengths and
      labeled non-speciation events.
    + *reconciled-tree_events.newick*: Rooted, reconciled tree with branch
      lengths in newick format. Tips are topiary UIDs. Internal nodes are labeled
      with evolutionary events.
    + *reconciliations/*: Directory with reconciliation information written out
      by GeneRax. See the `GeneRax documentation <generax-link_>`_ for details.
    + *reconciled-tree.newick*: Rooted, reconciled tree with branch lengths in
      newick format. Tips are topiary UIDs. Internal nodes are unlabeled.
    + *species-tree.newick*: Species tree downloaded from the Open Tree of Life.

  + Function:

    - `topiary.reconcile <topiary.generax.html#topiary.generax.reconcile.reconcile>`_
    - `topiary.df_to_species_tree <topiary.opentree.html#topiary.opentree.tree.df_to_species_tree>`_


+ :emph:`Infers ancestral sequences` on the reconciled tree using RAxML-NG and
  the maximum likelihood phylogenetic model. Gaps are inferred using parsimony
  via PastML. This tree an example of "summary-tree.pdf" generated by this
  step. The graph shows an example ancestor sequence summary. The text shows
  the corresponding reconstructed ancestor and its altAll versions as a fasta
  file.

  .. image:: _static/img/ali-to-anc/04_reconciled-tree-ancestors.svg
    :align: center
    :alt: ancestors drawn on reconciled gene/species tree
    :width: 75%

  :raw-html:`<br />`

  .. image:: _static/img/ali-to-anc/04_reconciled-tree-ancestors_plot.svg
    :align: center
    :alt: Graph showing ancestral quality
    :width: 85%

  :raw-html:`<br />`

  .. code-block::

    >anc12|avgPP:0.804|lnPP:-57.347|num_ambig:34|num_ambig_gaps:0
    MKVFFTLLFVLTLF-CS----G------GSKEWPTHTICNTSDLEVYYKSCDPL--Q-DVGLSISPCSKSMTENINIRVALILRQDIKELYLNLDLFIN----GLKVLSY--S--YPLCEPS-F-PKFTFCGRRKGELIYFEGPVKLGIQTIPQGEY-NVTLELFNEDN--Y--TIACVNITLI---------------S---R
    >anc12_altAll|avgPP:0.773|lnPP:-72.445|num_ambig:34|num_ambig_gaps:0
    MKIFFTLLLVLTLL-CT----G------ESKEWPTHTICNSSELEVYYRSCDPL--Q-DIGVSIEPCSKSLTENIDVRIALILRQDVKELYLDLDLYLN----GLKVLNY--S--YPLCEPS-F-PRFTFCGRKKGEMIYYEGPVKLSVQTLPKGEY-NVSLQLYNEDN--Y--TIACANVTLI---------------S---K

  + Output: :code:`04_reconciled-tree-ancestors/output/`

    + *summary-tree.pdf*: Drawing of reconciled tree with branch lengths,
      labeled non-speciation events, and ancestors with posterior probabilities
      as a color map.
    + *reconciled-tree_anc-pp.newick*: Rooted, reconciled tree with branch
      lengths in newick format. Tips are topiary UIDs. Internal nodes are
      ancestral posterior probabilities.
    + *reconciled-tree_anc-label.newick*: Rooted, reconciled tree with branch
      lengths in newick format. Tips are topiary UIDs. Internal nodes are
      ancestor names.
    + *reconciled-tree_ancestors/*: The *ancX.pdf* files are pdf summaries of
      ancestral reconstructions for each ancestor as labeled in *summary-tree.pdf*.
      An example is shown above. These graphs show the support for: the most
      likely reconstruction at each site (black circles), the next most likely
      reconstruction (red circles), the locations of gaps (gray shading), and
      the location of ambiguous gaps (purple dashes). Statistics at the top
      indicate ancestor quality.

  + Function: `topiary.generate_ancestors <topiary.raxml.html#topiary.raxml.ancestors.generate_ancestors>`_

+ :emph:`Generates bootstrap replicates from the ML gene tree`. This uses the
  ML gene tree (step *01_gene-tree*), :emph:`NOT` the reconciled tree as input and
  generates up to 1,000 bootstrap replicates. These will be fed into the next
  script, which calculates the final branch supports on the reconciled tree.

  .. image:: _static/img/ali-to-anc/05_gene-tree-bootstraps.svg
    :align: center
    :alt: bootstrap supports mapped to the ML gene tree
    :width: 75%

  + Output: :code:`05_gene-tree-bootstraps/output/`

    + *summary-tree.pdf*: Drawing of ML gene tree (:emph:`NOT` the reconciled
      tree) with branch lengths. Nodes are colored by bootstrap support. In
      figure, tree is rooted by midpoint.
    + *gene-tree_supports.newick*: ML gene tree with branch lengths in newick
      format. Tips are topiary UIDs. Internal nodes are labeled by bootstrap
      support.
    + *bootstrap_replicates/*: Directory has bootstrap replicate alignments and
      bootstrap trees that will be fed into the next step.

  + Function: `topiary.generate_bootstraps <file:///Users/harmsm/work/programming/git-clones/topiary/docs/build/html/topiary.raxml.html#topiary.raxml.bootstrap.generate_bootstraps>`_

.. _bootstrap-reconcile-pipeline:

Bootstrap Reconcile
===================

:emph:`Calculates branch supports for the reconciled tree`. Takes the bootstrap
replicate alignments and gene trees generated by RAxML-NG and feeds each one
into GeneRax for gene-species tree reconciliation. Uses the resulting set of
reconciled trees to calculate bootstrap branch supports on each of the nodes
on the reconciled tree.

  .. image:: _static/img/reconcile-bootstrap/06_reconciled-tree-bootstraps.svg
    :align: center
    :alt: final reconciled tree with events, ancestors, ancestor supports, and branch supports
    :width: 75%

  + Output: :code:`06_reconciled-tree-bootstraps/output/`

    + *summary-tree.pdf*: Drawing of rooted reconciled tree with with branch
      lengths. Nodes are labeled with non-speciation evolutionary events,
      ancestor names, ancestor posterior probabilities (as color map), and
      branch supports (as color map).
    + *bootstrap-convergence-report.csv*: Spreadsheet showing results of a
      convergence test for the branch supports on the bootstrap tree. The key
      columns are the number of trees included (:code:`trees`) and the number of
      tree permutations converged (:code:`perms_below_cutoff`). See
      `Pattengale et. al. <http://www.liebertpub.com/doi/10.1089/cmb.2009.0179>`_ for details.
    + *reconciled-tree_supports.newick*: Rooted, reconciled tree with branch
      lengths in newick format. Tips are topiary UIDs. Internal nodes branch

  + `topiary.reconcile <topiary.generax.html#topiary.generax.reconcile.reconcile>`_.
    This is called with the :code:`bootstrap = True` flag.
  