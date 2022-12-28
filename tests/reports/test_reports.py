
from topiary.reports.reports import _find_directories
from topiary.reports.reports import tree_report
from topiary.reports.reports import pipeline_report

import os
import shutil

def test__find_directories(tmpdir,small_phylo):

    cwd = os.getcwd()
    os.chdir(tmpdir)

    # Complete directory
    os.mkdir("test1")
    dirs_to_copy = ["00_find-model",
                    "01_gene-tree",
                    "02_gene-tree-ancestors",
                    "03_reconciled-tree",
                    "04_reconciled-tree-ancestors",
                    "05_gene-tree-bootstraps",
                    "06_reconciled-tree-bootstraps"]

    for d in dirs_to_copy:
        shutil.copytree(small_phylo[d],os.path.join("test1",d))
    
    out = _find_directories("test1")

    assert out['model'] == 'test1/00_find-model'
    assert out['gene']['anc'] == 'test1/02_gene-tree-ancestors'
    assert out['gene']['tree'] == 'test1/05_gene-tree-bootstraps'
    assert out['reconciled']['anc'] == 'test1/04_reconciled-tree-ancestors'
    assert out['reconciled']['tree'] == 'test1/06_reconciled-tree-bootstraps'

    # -------------------------------------------------------------------------
    # No reconciled bootstraps

    os.mkdir("test2")
    dirs_to_copy = ["00_find-model",
                    "01_gene-tree",
                    "02_gene-tree-ancestors",
                    "03_reconciled-tree",
                    "04_reconciled-tree-ancestors",
                    "05_gene-tree-bootstraps"]

    for d in dirs_to_copy:
        shutil.copytree(small_phylo[d],os.path.join("test2",d))
    
    out = _find_directories("test2")

    assert out['model'] == 'test2/00_find-model'
    assert out['gene']['anc'] == 'test2/02_gene-tree-ancestors'
    assert out['gene']['tree'] == 'test2/05_gene-tree-bootstraps'
    assert out['reconciled']['anc'] == 'test2/04_reconciled-tree-ancestors'
    assert out['reconciled']['tree'] == 'test2/04_reconciled-tree-ancestors'

    # -------------------------------------------------------------------------
    # No gene tree bootstraps

    os.mkdir("test3")
    dirs_to_copy = ["00_find-model",
                    "01_gene-tree",
                    "02_gene-tree-ancestors",
                    "03_reconciled-tree",
                    "04_reconciled-tree-ancestors"]

    for d in dirs_to_copy:
        shutil.copytree(small_phylo[d],os.path.join("test3",d))
    
    out = _find_directories("test3")

    assert out['model'] == 'test3/00_find-model'
    assert out['gene']['anc'] == 'test3/02_gene-tree-ancestors'
    assert out['gene']['tree'] == 'test3/02_gene-tree-ancestors'
    assert out['reconciled']['anc'] == 'test3/04_reconciled-tree-ancestors'
    assert out['reconciled']['tree'] == 'test3/04_reconciled-tree-ancestors'

    # -------------------------------------------------------------------------
    # No reconciled ancestors
    
    os.mkdir("test4")
    dirs_to_copy = ["00_find-model",
                    "01_gene-tree",
                    "02_gene-tree-ancestors",
                    "03_reconciled-tree"]

    for d in dirs_to_copy:
        shutil.copytree(small_phylo[d],os.path.join("test4",d))
    
    out = _find_directories("test4")

    assert out['model'] == 'test4/00_find-model'
    assert out['gene']['anc'] == 'test4/02_gene-tree-ancestors'
    assert out['gene']['tree'] == 'test4/02_gene-tree-ancestors'
    assert out['reconciled']['anc'] is None
    assert out['reconciled']['tree'] == 'test4/03_reconciled-tree'


    # -------------------------------------------------------------------------
    # No reconciled tree
    
    os.mkdir("test5")
    dirs_to_copy = ["00_find-model",
                    "01_gene-tree",
                    "02_gene-tree-ancestors"]

    for d in dirs_to_copy:
        shutil.copytree(small_phylo[d],os.path.join("test5",d))
    
    out = _find_directories("test5")

    assert out['model'] == 'test5/00_find-model'
    assert out['gene']['anc'] == 'test5/02_gene-tree-ancestors'
    assert out['gene']['tree'] == 'test5/02_gene-tree-ancestors'
    assert out['reconciled']['anc'] is None
    assert out['reconciled']['tree'] is None

    # -------------------------------------------------------------------------
    # No gene tree ancestors
    
    os.mkdir("test6")
    dirs_to_copy = ["00_find-model",
                    "01_gene-tree"]

    for d in dirs_to_copy:
        shutil.copytree(small_phylo[d],os.path.join("test6",d))
    
    out = _find_directories("test6")

    assert out['model'] == 'test6/00_find-model'
    assert out['gene']['anc'] is None
    assert out['gene']['tree'] == 'test6/01_gene-tree'
    assert out['reconciled']['anc'] is None
    assert out['reconciled']['tree'] is None

    # -------------------------------------------------------------------------
    # No gene tree 
    
    os.mkdir("test7")
    dirs_to_copy = ["00_find-model"]

    for d in dirs_to_copy:
        shutil.copytree(small_phylo[d],os.path.join("test7",d))
    
    out = _find_directories("test7")

    assert out['model'] == 'test7/00_find-model'
    assert out['gene']['anc'] is None
    assert out['gene']['tree'] is None
    assert out['reconciled']['anc'] is None
    assert out['reconciled']['tree'] is None

    # -------------------------------------------------------------------------
    # Nothing at all
    
    os.mkdir("test8")
    dirs_to_copy = []

    for d in dirs_to_copy:
        shutil.copytree(small_phylo[d],os.path.join("test8",d))
    
    out = _find_directories("test8")

    assert out['model'] is None
    assert out['gene']['anc'] is None
    assert out['gene']['tree'] is None
    assert out['reconciled']['anc'] is None
    assert out['reconciled']['tree'] is None

    # -------------------------------------------------------------------------
    # Single calc dir, gene tree

    out = _find_directories(small_phylo["01_gene-tree"])

    assert out['model'] is None
    assert out['gene']['anc'] is None
    assert os.path.split(out['gene']['tree'])[-1] == "01_gene-tree"
    assert out['reconciled']['anc'] is None
    assert out['reconciled']['tree'] is None

    # -------------------------------------------------------------------------
    # Single calc dir, model

    out = _find_directories(small_phylo["00_find-model"])

    assert os.path.split(out['model'])[-1] == "00_find-model"
    assert out['gene']['anc'] is None
    assert out["gene"]["tree"] is None
    assert out['reconciled']['anc'] is None
    assert out['reconciled']['tree'] is None

    os.chdir(cwd)


    

def test_tree_report():
    pass 

def test_pipeline_report():
    pass