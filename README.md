# topiary

Python framework for doing ancestral sequence reconstruction and other
phylogenetics using pandas dataframes and ete3 trees as the primary data
structures.

*Under heavy development.*

### Current installation instructions:

On a command line/conda terminal:

```
git clone https://github.com/harmslab/topiary.git
cd topiary
python setup.py install

conda install -c bioconda muscle
conda install -c bioconda generax
conda install -c bioconda raxml-ng
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

If the BLAST bit of the pipeline doesn't work, you may need to install a BLAST binary.
You can download these from the [NCBI](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
