# topiary

Python framework for doing ancestral sequence reconstruction using pandas dataframes and ete3 trees as the primary data structures.

*Under heavy development.*

### Current installation instructions:

On a command line/conda terminal:

```
conda install -c bioconda muscle
conda install -c bioconda generax
conda install -c bioconda raxml-ng
conda install -c bioconda blast

git clone https://github.com/harmslab/topiary.git
cd topiary
python setup.py install
```

To see an example that uses the code, use the terminal to go into a new
directory (not the one with topiary) and run:

```
git clone https://github.com/harmslab/asr-protocol.git
cd asr-protocol/asr-protocol
jupyter-lab
```

Open the `quickstart.ipynb` notebook. You should be able to run an ancestral
sequence reconstruction calculation on a test dataset by running through this
notebook.


### Pipeline and citations:

+ Muscle [Edgar RC (2021) *bioRxiv* https://doi.org/10.1101/2021.06.20.449169](https://doi.org/10.1101/2021.06.20.449169)
+ Tree and ancestor reconstruction using RaxML-NG [Kozlov et al (2019) *Bioinformatics* https://doi.org/10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)
+ Gene/species tree reconciliation with [GeneRax](https://github.com/BenoitMorel/GeneRax)
[Morel et al (2020) *MBE* https://doi.org/10.1093/molbev/msaa141](https://doi.org/10.1093/molbev/msaa141)
+ Ancestral gap assignment using parsimony, as implemented in PastML 
[Ishikawa et al (2019) *MBE* https://doi.org/10.1093/molbev/msz131](https://doi.org/10.1093/molbev/msz131)
+ Species tree access using [https://tree.opentreeoflife.org/](https://tree.opentreeoflife.org/).




