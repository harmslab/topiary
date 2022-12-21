
.. role:: emph

.. role:: raw-html(raw)
    :format: html

.. include:: links.rst

.. _taxonomic-scope-doc:

===========================
Determining taxonomic scope
===========================

The following steps can be used to determine the taxonomic scope for a protein 
of interest. We recommend doing this for non-microbial proteins. For microbial
proteins, the species tree is poorly defined; therefore, topiary automatically
sets the scope to microbial proteins. 

#. BLAST known protein sequences against the 
   `non-redundant clustered <nr-clustered-database_>`_ BLAST database, setting
   the :code:`Max target sequences` parameter to 1,000 or more.

   .. image:: _static/img/find-taxonomic-scope/step-1.png
     :align: center
     :alt: BLAST sequence against non-redundant clustered
     :height: 300

   :raw-html:`<br />`
#. Take a few representative sequences from the most divergent species in the
   outputs. These can be selected using the "Taxonomy" tab.

   .. image:: _static/img/find-taxonomic-scope/step-2.png
     :align: center
     :alt: BLAST sequence against non-redundant clustered
     :height: 300

   :raw-html:`<br />`
#. BLAST the divergent sequences back against the NCBI 
   `non-redundant <nr-database_>`_ database, limiting the search to the
   species from which you took your known sequences. (If, for example, you
   used a human sequence as your starting point, you would limit this
   "reciprocal" query to  *Homo sapiens*.)
  
   .. image:: _static/img/find-taxonomic-scope/step-3.png
     :align: center
     :alt: BLAST sequence against non-redundant clustered
     :height: 300
  
   :raw-html:`<br />`
#. If this BLAST search pulls up your starting protein as a top hit, it is
   good evidence that the species from which the sequence came has the protein
   and should be included in the taxonomic scope. 
 
   .. image:: _static/img/find-taxonomic-scope/step-4.png
     :align: center
     :alt: BLAST sequence against non-redundant clustered
     :height: 300
 
   :raw-html:`<br />`