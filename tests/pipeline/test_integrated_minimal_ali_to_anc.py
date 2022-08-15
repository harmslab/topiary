import pytest

import topiary

import os

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_integrated_minimal_ali_to_anc(tiny_phylo,tmpdir):
    """
    Full test of ali_to_anc pipeline without any complexity of arg checking
    etc. Goal is to catch major problems or changes to core functionality
    linking pieces together. The generax steps are all set to have thread = 1
    to avoid any mpi complexity. We also pass in a known species tree to
    avoid any interaction with open tree of life. 
    """

    df = tiny_phylo["initial-input/dataframe.csv"]
    species_tree = tiny_phylo["initial-input/species-tree.newick"]

    def _check_out_files(tiny_phylo,out_dir):

        key = "/".join([out_dir,"output","*"])
        expected_files = tiny_phylo[key]
        for e in expected_files:
            assert os.path.exists(os.path.join(out_dir,"output",
                                               os.path.basename(e)))


    current_dir = os.getcwd()
    os.chdir(tmpdir)

    # -------------------------------------------------------------------------
    # find best model

    topiary.find_best_model(df=df,
                            calc_dir="00_find-best-model",
                            seed=12345,
                            model_matrices=["LG","JTT"],
                            model_rates=None,
                            model_freqs=None,
                            model_invariant=None)

    _check_out_files(tiny_phylo,"00_find-best-model")

    # -------------------------------------------------------------------------
    # generate ml tree

    topiary.generate_ml_tree(previous_dir="00_find-best-model",
                             calc_dir="01_gene-tree")

    _check_out_files(tiny_phylo,"01_gene-tree")

    # -------------------------------------------------------------------------
    # reconcile

    topiary.reconcile(previous_dir="01_gene-tree/",
                      calc_dir="02_reconcile",
                      species_tree=species_tree,
                      num_threads=1)

    _check_out_files(tiny_phylo,"02_reconcile")

    # -------------------------------------------------------------------------
    # infer ancestors

    topiary.generate_ancestors(previous_dir="02_reconcile/",
                               calc_dir="03_ancestors")

    _check_out_files(tiny_phylo,"03_ancestors")

    # -------------------------------------------------------------------------
    # gene tree bootstraps

    topiary.generate_bootstraps(previous_dir="03_ancestors/",
                                calc_dir="04_bootstraps")

    _check_out_files(tiny_phylo,"04_bootstraps")

    # -------------------------------------------------------------------------
    # reconcile bootstraps

    topiary.reconcile(previous_dir="04_bootstraps/",
                      calc_dir="05_reconcile-bootstraps",
                      num_threads=1,
                      bootstrap=True)

    _check_out_files(tiny_phylo,"05_reconcile-bootstraps")

    os.chdir(current_dir)
