# topiary: python framework for doing ancestral sequence reconstruction

![Testing status](https://github.com/harmslab/topiary/actions/workflows/python-app.yml/badge.svg) ![Coverage](docs/badges/coverage-badge.svg) ![Number of tests](docs/badges/tests-badge.svg) ![Documentation Status](https://readthedocs.org/projects/topiary-asr/badge/?version=latest)

### [Installation](https://topiary-asr.readthedocs.io/en/latest/installation.html)

### [Documentation](https://topiary-asr.readthedocs.io/en/latest/)

### Try it out on Google Colab

+ [Go from a few initial sequences to a full alignment](https://githubtocolab.com/harmslab/topiary-examples/blob/main/notebooks/seed-to-alignment.ipynb)
+ [Build a phylogenetic tree and reconstruct ancestors](https://githubtocolab.com/harmslab/topiary-examples/blob/main/notebooks/alignment-to-ancestors.ipynb)

![ASR pipeline](docs/source/_static/img/asr-pipeline-01.png)

#### Features

+ *Automatic.* Performs sequence database construction, quality
  control, multiple sequence alignment, tree construction, gene/species tree
  reconciliation, and ancestral reconstruction with minimal user input.
+ *Species aware.* Integrates with the [Open Tree of Life](https://tree.opentreeoflife.org/opentree/argus/opentree13.4)
  database, improving selection of sequences and tree/ancestor inference.
+ *Human-oriented.* Users prepare input as spreadsheets, not
  arcane text files. Outputs are spreadsheets, clean fasta files, pdf trees,
  and graphical summaries of ancestor quality.
+ *Flexible.* Use as a command line program or do custom analyses
  and plotting with the API.
+ *Modern.* Topiary is built around a collection of modern,
  actively-supported, phylogenetic software tools:

  + [OpenTree](https://opentree.readthedocs.io/en/latest/)
  + [muscle 5](https://www.drive5.com/muscle/)
  + [RAxML-NG](https://github.com/amkozlov/raxml-ng)
  + [GeneRax](https://github.com/BenoitMorel/GeneRax)
  + [PastML](https://pastml.pasteur.fr)
  + [toytree](https://toyplot.readthedocs.io/)

#### Citing topiary
Orlandi KN, Phillips SR, Sailer ZR, Harman JL, Harms MJ. "Topiary: pruning the
manual labor from ancestral sequence reconstruction" (2022) *Protein Science* 
[doi: 10.1002/pro.4551](https://onlinelibrary.wiley.com/doi/10.1002/pro.4551)
