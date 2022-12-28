import pytest
import topiary

from topiary.io.tree import read_tree
from topiary.io.tree import _map_tree_to_tree
from topiary.io.tree import load_trees
from topiary.io.tree import write_trees

import os
import re

def test_read_tree():

    pass

def test__map_tree_to_tree():

    pass

def test_load_trees():

    pass

def test_write_trees(small_phylo,tmpdir):

    cwd = os.getcwd()
    os.chdir(tmpdir)

    # Reconciled trees
    T = load_trees(small_phylo["06_reconciled-tree-bootstraps/output"],prefix="reconciled")
    newick = write_trees(T)
    split_trees = newick.split(";")[:-1]
    assert len(split_trees) == 4
    
    expected = ["100","S","a5","0.999999"]
    for i in range(4):
        last_value = split_trees[i].split(")")[-2].split(":")[0]
        assert last_value == expected[i]
        
    # gene trees
    T = load_trees(small_phylo["06_reconciled-tree-bootstraps/output"],prefix="gene")
    newick = write_trees(T)
    split_trees = newick.split(";")[:-1]
    assert len(split_trees) == 3
    
    expected = ["100","a10","0.997959"]
    for i in range(3):
        last_value = split_trees[i].split(")")[-2].split(":")[0]
        assert last_value == expected[i]
        
    # Only anc_pp
    T = load_trees(small_phylo["06_reconciled-tree-bootstraps/output"],prefix="reconciled")
    newick = write_trees(T,
                         anc_pp=True,
                         anc_label=False,
                         bs_support=False,
                         event=False)
    split_trees = newick.split(";")[:-1]
    assert len(split_trees) == 1
    last_value = split_trees[0].split(")")[-2].split(":")[0]
    assert last_value == "0.999999"

    # Only anc_label
    T = load_trees(small_phylo["06_reconciled-tree-bootstraps/output"],prefix="reconciled")
    newick = write_trees(T,
                         anc_pp=False,
                         anc_label=True,
                         bs_support=False,
                         event=False)
    split_trees = newick.split(";")[:-1]
    assert len(split_trees) == 1
    last_value = split_trees[0].split(")")[-2].split(":")[0]
    assert last_value == "a5"

    # Only bs_support
    T = load_trees(small_phylo["06_reconciled-tree-bootstraps/output"],prefix="reconciled")
    newick = write_trees(T,
                         anc_pp=False,
                         anc_label=False,
                         bs_support=True,
                         event=False)
    split_trees = newick.split(";")[:-1]
    assert len(split_trees) == 1
    last_value = split_trees[0].split(")")[-2].split(":")[0]
    assert last_value == "100"

    # Only event
    T = load_trees(small_phylo["06_reconciled-tree-bootstraps/output"],prefix="reconciled")
    newick = write_trees(T,
                         anc_pp=False,
                         anc_label=False,
                         bs_support=False,
                         event=True)
    split_trees = newick.split(";")[:-1]
    assert len(split_trees) == 1
    last_value = split_trees[0].split(")")[-2].split(":")[0]
    assert last_value == "S"

    # Write to file
    T = load_trees(small_phylo["06_reconciled-tree-bootstraps/output"],prefix="reconciled")
    newick = write_trees(T,out_file="test.newick")
    f = open("test.newick","r")
    lines = f.readlines()
    f.close()
    assert len(lines) == 4
    assert sum([line.strip()[-1] == ";" for line in lines]) == 4

    # Try to write to existing
    with pytest.raises(FileExistsError):
        newick = write_trees(T,out_file="test.newick")
                
    # Overwrite
    newick = write_trees(T,out_file="test.newick",overwrite=True)

    # Try to write onto a directory
    os.mkdir("test")
    with pytest.raises(FileExistsError):
        newick = write_trees(T,out_file="test",overwrite=True)

    # name_dict
    T = load_trees(small_phylo["06_reconciled-tree-bootstraps/output"],prefix="reconciled")

    # Build simple name_dict
    name_dict = {}
    for n in T.traverse():
        if n.is_leaf():
            name_dict[n.name] = re.sub("x","Y",n.name)
    
    # Make sure replacement happened
    newick = write_trees(T,name_dict=name_dict)
    found = re.search("Y",newick)
    assert found is not None
    found = re.search("x",newick)
    assert found is None
    
    os.chdir(cwd)