import pytest

import topiary

import os

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_integrated_minimal_ali_to_anc(tiny_phylo,tmpdir):
    """
    Full test of ali_to_anc pipeline without any complexity of arg checking
    etc. Goal is to catch major problems or changes to core functionality
    linking pieces together. We pass in a known species tree to
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

    topiary.generate_ml_tree(prev_calculation="00_find-best-model",
                             calc_dir="01_gene-tree")

    _check_out_files(tiny_phylo,"01_gene-tree")

    # -------------------------------------------------------------------------
    # reconcile
    # if num_threads > 1 can run into errors because this calculation is so
    # small. generax crashes because it tries to read file on one thread when
    # it hasn't yet finished writing on another thread. not really thread safe
    # :grimace:
    topiary.reconcile(prev_calculation="01_gene-tree",
                      calc_dir="02_reconcile",
                      species_tree=species_tree,
                      num_threads=1)

    _check_out_files(tiny_phylo,"02_reconcile")

    # -------------------------------------------------------------------------
    # infer ancestors

    topiary.generate_ancestors(prev_calculation="02_reconcile",
                               calc_dir="03_ancestors")

    _check_out_files(tiny_phylo,"03_ancestors")

    # -------------------------------------------------------------------------
    # gene tree bootstraps

    topiary.generate_bootstraps(prev_calculation="03_ancestors",
                                calc_dir="04_bootstraps",
                                num_threads=-1)

    _check_out_files(tiny_phylo,"04_bootstraps")

    # -------------------------------------------------------------------------
    # reconcile bootstraps

    topiary.reconcile(prev_calculation=tiny_phylo["04_bootstraps_toy"],
                      calc_dir="05_reconcile-bootstraps",
                      bootstrap=True)

    # Because we are using the topy bootstraps rather than real bootstrap from
    # last step, only check for primary expected output; don't worry about
    # other files that could be in that directory.
    assert os.path.isfile(os.path.join("05_reconcile-bootstraps",
                                       "output",
                                       "reconciled-tree_supports.newick"))


    os.chdir(current_dir)
