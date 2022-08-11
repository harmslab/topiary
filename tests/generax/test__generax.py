import pytest

import topiary
from topiary.generax._generax import _get_link_dict
from topiary.generax._generax import setup_generax
from topiary.generax._generax import run_generax
from topiary.generax._generax import GENERAX_BINARY

import ete3

import numpy as np
import pandas as pd

import os
import json
import shutil
import pathlib

def test__get_link_dict():

    df = pd.DataFrame({"ott":["A","B","C","D"],
                       "uid":["0","1","2","3"]})
    T = ete3.Tree("((0,1),(2,3));")

    link_dict, uid_in_gene_tree = _get_link_dict(df,T)
    assert np.array_equal(link_dict["A"],["0"])
    assert np.array_equal(link_dict["B"],["1"])
    assert np.array_equal(link_dict["C"],["2"])
    assert np.array_equal(link_dict["D"],["3"])

    uid_in_gene_tree.sort()
    assert np.array_equal(uid_in_gene_tree,["0","1","2","3"])


    # degeneracy

    df = pd.DataFrame({"ott":["A","A","B","B"],
                       "uid":["0","1","2","3"]})
    T = ete3.Tree("((0,1),(2,3));")

    link_dict, uid_in_gene_tree = _get_link_dict(df,T)
    A = link_dict["A"]
    A.sort()
    assert np.array_equal(A,["0","1"])

    B = link_dict["B"]
    B.sort()
    assert np.array_equal(B,["2","3"])

    uid_in_gene_tree.sort()
    assert np.array_equal(uid_in_gene_tree,["0","1","2","3"])

    # only one
    df = pd.DataFrame({"ott":["A","A","A","A"],
                       "uid":["0","1","2","3"]})
    T = ete3.Tree("((0,1),(2,3));")

    link_dict, uid_in_gene_tree = _get_link_dict(df,T)
    A = link_dict["A"]
    A.sort()
    assert np.array_equal(A,["0","1","2","3"])
    with pytest.raises(KeyError):
        link_dict["B"]

    uid_in_gene_tree.sort()
    assert np.array_equal(uid_in_gene_tree,["0","1","2","3"])

    # missing uid from tree
    df = pd.DataFrame({"ott":["A","A","B","B"],
                       "uid":["0","1","2","3"]})
    T = ete3.Tree("((0,1),(2));")

    link_dict, uid_in_gene_tree = _get_link_dict(df,T)
    A = link_dict["A"]
    A.sort()
    assert np.array_equal(A,["0","1"])

    B = link_dict["B"]
    B.sort()
    assert np.array_equal(B,["2"])

    uid_in_gene_tree.sort()
    assert np.array_equal(uid_in_gene_tree,["0","1","2"])



def test_setup_generax(generax_data,tmpdir):

    # Make sure it reads in previous run without error
    prev_ml_dir = generax_data["prev-ml-run"]
    ml_only = os.path.join(prev_ml_dir,"01_ml-tree")
    ml_bs = os.path.join(prev_ml_dir,"04_bootstraps")

    for i, d in enumerate([ml_only,ml_bs]):

        out = os.path.join(d,"output")

        df = topiary.read_dataframe(os.path.join(out,"dataframe.csv"))
        gene_tree = os.path.join(out,"tree.newick")

        f = open(os.path.join(out,"run_parameters.json"),'r')
        params = json.load(f)
        f.close()
        model = params['model']

        out_dir = os.path.join(tmpdir,f"out_{i}")

        keep_mask = setup_generax(df=df,
                                  gene_tree=gene_tree,
                                  model=model,
                                  out_dir=out_dir)

    # -----------------------------------------------------------------------
    # Check a toy file with expected outputs

    toy_inputs = os.path.join(generax_data["toy-input"],"toy-ml","output")
    toy_out = os.path.join(tmpdir,"toy-ml-setup-generax")
    df = topiary.read_dataframe(os.path.join(toy_inputs,"dataframe.csv"))
    gene_tree = os.path.join(toy_inputs,"tree.newick")

    f = open(os.path.join(toy_inputs,"run_parameters.json"),'r')
    params = json.load(f)
    f.close()
    model = params['model']

    keep_mask = setup_generax(df=df,
                              gene_tree=gene_tree,
                              model=model,
                              out_dir=toy_out)

    assert np.array_equal(keep_mask,np.ones(len(keep_mask)))

    # make sure species tree has expecetd tips
    st_out = os.path.join(toy_out,"species_tree.newick")
    assert os.path.isfile(st_out)
    T = ete3.Tree(st_out)
    leaf_names = set(T.get_leaf_names())
    expected = set(["ott276534","ott913930","ott542509","ott271571"])
    assert len(leaf_names - expected) == 0

    # check link file
    def _read_link(link_file):

        out_dict = {}
        with open(link_file) as f:
            for line in f:
                key = line.split(":")[0].strip()
                values = [v.strip() for v in line.split(":")[1].split(";")]
                values.sort()

                out_dict[key] = tuple(values)

        return out_dict

    out_link = _read_link(os.path.join(toy_out,"mapping.link"))
    expected_link = _read_link(os.path.join(toy_inputs,"..","expected-output","mapping.link"))
    for k in out_link:
        assert out_link[k] == expected_link[k]

    # Make sure the alignments are equivalent
    f = open(os.path.join(toy_out,"alignment.phy"))
    out_ali = [line.strip() for line in f.readlines()]
    f.close()

    f = open(os.path.join(toy_inputs,"..","expected-output","alignment.phy"))
    expected_ali = [line.strip() for line in f.readlines()]
    f.close()

    np.array_equal(out_ali,expected_ali)

    # Make sure control files are same
    f = open(os.path.join(toy_out,"control.txt"))
    out_ctrl = [line.strip() for line in f.readlines()]
    f.close()

    f = open(os.path.join(toy_inputs,"..","expected-output","control.txt"))
    expected_ctrl= [line.strip() for line in f.readlines()]
    f.close()

    np.array_equal(out_ctrl,expected_ctrl)

    # Make sure control files are same
    f = open(os.path.join(toy_out,"gene_tree.newick"))
    out_gt = [line.strip() for line in f.readlines()]
    f.close()

    f = open(os.path.join(toy_inputs,"..","expected-output","gene_tree.newick"))
    expected_gt= [line.strip() for line in f.readlines()]
    f.close()

    np.array_equal(out_gt,expected_gt)

    # ---
    # check pass-in capabilities

    # keep_mask alone will get wiped out (noted on docs)
    shutil.rmtree(toy_out)

    in_keep_mask = np.array([0,1,0,1,0,1,0,1],dtype=bool)
    keep_mask = setup_generax(df=df,
                              gene_tree=gene_tree,
                              model=model,
                              out_dir=toy_out,
                              keep_mask=in_keep_mask)

    assert np.array_equal(np.ones(len(keep_mask)),keep_mask)


    shutil.rmtree(toy_out)

    # Send in wacky keep mask and gene_tree as newick. (That file is just to
    # make sure it's actually being copied in; contents don't really matter for
    # this test; this function does not validate file type/content.)
    in_keep_mask = np.array([0,1,0,1,0,1,0,1],dtype=bool)
    in_mapping_link = os.path.join(toy_inputs,"..","expected-output","gene_tree.newick")
    keep_mask = setup_generax(df=df,
                              gene_tree=gene_tree,
                              model=model,
                              out_dir=toy_out,
                              mapping_link_file=in_mapping_link,
                              keep_mask=in_keep_mask)

    assert np.array_equal(in_keep_mask,keep_mask)

    # Make sure  files are same
    f = open(os.path.join(toy_out,"mapping.link"))
    out_link = [line.strip() for line in f.readlines()]
    f.close()

    f = open(os.path.join(in_mapping_link))
    expected_link = [line.strip() for line in f.readlines()]
    f.close()

    np.array_equal(out_link,expected_link)


    shutil.rmtree(toy_out)

    # Send in file for species_tree_file too
    in_keep_mask = np.array([0,1,0,1,0,1,0,1],dtype=bool)
    in_mapping_link = os.path.join(toy_inputs,"..","expected-output","gene_tree.newick")
    keep_mask = setup_generax(df=df,
                              gene_tree=gene_tree,
                              model=model,
                              out_dir=toy_out,
                              species_tree_file=in_mapping_link,
                              mapping_link_file=in_mapping_link,
                              keep_mask=in_keep_mask)

    assert np.array_equal(in_keep_mask,keep_mask)

    # Make sure files are same
    f = open(os.path.join(toy_out,"mapping.link"))
    out_link = [line.strip() for line in f.readlines()]
    f.close()

    f = open(os.path.join(in_mapping_link))
    expected_link = [line.strip() for line in f.readlines()]
    f.close()

    np.array_equal(out_link,expected_link)

    f = open(os.path.join(toy_out,"species_tree.newick"))
    out_st = [line.strip() for line in f.readlines()]
    f.close()

    np.array_equal(out_st,expected_link)

    shutil.rmtree(toy_out)

    # Send in file for species_tree_file too
    in_keep_mask = np.array([0,1,0,1,0,1,0,1],dtype=bool)
    in_mapping_link = os.path.join(toy_inputs,"..","expected-output","gene_tree.newick")
    keep_mask = setup_generax(df=df,
                              gene_tree=gene_tree,
                              model=model,
                              out_dir=toy_out,
                              control_file=in_mapping_link,
                              species_tree_file=in_mapping_link,
                              mapping_link_file=in_mapping_link,
                              keep_mask=in_keep_mask)

    assert np.array_equal(in_keep_mask,keep_mask)

    # Make sure files are same
    f = open(os.path.join(toy_out,"mapping.link"))
    out_link = [line.strip() for line in f.readlines()]
    f.close()

    f = open(os.path.join(in_mapping_link))
    expected_link = [line.strip() for line in f.readlines()]
    f.close()

    np.array_equal(out_link,expected_link)

    f = open(os.path.join(toy_out,"species_tree.newick"))
    out_st = [line.strip() for line in f.readlines()]
    f.close()

    np.array_equal(out_st,expected_link)

    f = open(os.path.join(toy_out,"control.txt"))
    out_ctrl = [line.strip() for line in f.readlines()]
    f.close()

    np.array_equal(out_ctrl,expected_link)

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_run_generax(generax_data,tmpdir):

    # Validate can run from both ml and ml_bootstraps output
    # Make sure it reads in previous run without error
    prev_ml_dir = generax_data["prev-ml-run"]
    ml_only = os.path.join(prev_ml_dir,"01_ml-tree")
    ml_bs = os.path.join(prev_ml_dir,"04_bootstraps")

    for i, d in enumerate([ml_only,ml_bs]):

        out = os.path.join(d,"output")

        df = topiary.read_dataframe(os.path.join(out,"dataframe.csv"))
        gene_tree = os.path.join(out,"tree.newick")

        f = open(os.path.join(out,"run_parameters.json"),'r')
        params = json.load(f)
        f.close()
        model = params['model']

        out_dir = os.path.join(tmpdir,f"out_{i}")
        if os.path.isdir(out_dir):
            shutil.rmtree(out_dir)

        keep_mask = setup_generax(df=df,
                                  gene_tree=gene_tree,
                                  model=model,
                                  out_dir=out_dir)

        cmd = run_generax(run_directory=out_dir,
                          allow_horizontal_transfer=True,
                          num_threads=1,
                          generax_binary=GENERAX_BINARY,
                          write_to_script="run_generax.sh")

        # Make sure command is being constructed correctly and that it is being
        # written properly to the script
        assert cmd == "generax --families control.txt --species-tree species_tree.newick --prefix result --rec-model UndatedDTL"
        f = open(os.path.join(out_dir,"run_generax.sh"))
        content = f.read()
        f.close()
        assert cmd == content.strip()

        # Get horizontal transfer diff
        cmd = run_generax(run_directory=out_dir,
                          allow_horizontal_transfer=False,
                          num_threads=1,
                          generax_binary=GENERAX_BINARY,
                          write_to_script="run_generax.sh")

        assert cmd == "generax --families control.txt --species-tree species_tree.newick --prefix result --rec-model UndatedDL"

        # Validate num_threads
        cmd = run_generax(run_directory=out_dir,
                          allow_horizontal_transfer=True,
                          num_threads=2,
                          generax_binary=GENERAX_BINARY,
                          write_to_script="run_generax.sh")

        assert cmd == "mpirun -np 2 generax --families control.txt --species-tree species_tree.newick --prefix result --rec-model UndatedDTL"

        with pytest.raises(ValueError):
            cmd = run_generax(run_directory=out_dir,
                              allow_horizontal_transfer=True,
                              num_threads=1000000,
                              generax_binary=GENERAX_BINARY,
                              write_to_script="run_generax.sh")

        with pytest.raises(ValueError):
            cmd = run_generax(run_directory=out_dir,
                              allow_horizontal_transfer=True,
                              num_threads=1,
                              generax_binary="not_a_binary",
                              write_to_script="run_generax.sh")

        with pytest.raises(FileNotFoundError):
            cmd = run_generax(run_directory="not_an_out_dir",
                              allow_horizontal_transfer=True,
                              num_threads=1,
                              generax_binary=GENERAX_BINARY,
                              write_to_script="run_generax.sh")

        # Send in file for out dir
        test_file = os.path.join(tmpdir,"test_file")
        pathlib.Path(test_file).touch()

        with pytest.raises(ValueError):
            cmd = run_generax(run_directory=test_file,
                              allow_horizontal_transfer=True,
                              num_threads=1,
                              generax_binary=GENERAX_BINARY,
                              write_to_script="run_generax.sh")

    # -------------------------------------------------------------------------
    # Validate one reconcilation toy data. In this dataset, the gene tree has
    # a (human,mouse),(lemur,rat) topology for one of the paralogs. In test,
    # make sure the correct topology and evolutionary events are recovered.

    # Load in toy data

    input_dir = os.path.join(generax_data["toy-input"],"toy-ml","output")

    df = topiary.read_dataframe(os.path.join(input_dir,"dataframe.csv"))
    gene_tree = os.path.join(input_dir,"tree_wrong.newick")

    f = open(os.path.join(input_dir,"run_parameters.json"),'r')
    params = json.load(f)
    f.close()
    model = params['model']

    species_tree_file = os.path.join(input_dir,"species_tree.newick")

    run_directory = os.path.join(tmpdir,"real-generax-toy-test")
    if os.path.exists(run_directory):
        shutil.rmtree(run_directory)

    # Make directory in which to do the calculation
    keep_mask = setup_generax(df=df,
                              gene_tree=gene_tree,
                              model=model,
                              species_tree_file=species_tree_file,
                              out_dir=run_directory)

    # Run the calculation
    cmd = run_generax(run_directory=run_directory,
                      allow_horizontal_transfer=True,
                      num_threads=1,
                      generax_binary=GENERAX_BINARY)

    result_dir = os.path.join(run_directory,"result","results","reconcile")

    assert os.path.isfile(os.path.join(result_dir,"geneTree.newick"))

    # Make sure the reconciliation events tree exists and that it properly finds
    # a single duplication, six speciations, and no other evolutionary events
    reconcile_tree = os.path.join(run_directory,"result",
                                  "reconciliations","reconcile_events.newick")
    T = ete3.Tree(reconcile_tree,format=1)

    events = {}
    for n in T.traverse():
        if not n.is_leaf():
            try:
                events[n.name] += 1
            except KeyError:
                events[n.name] = 1

    # Assert correct topology is inferred
    assert len(T.get_common_ancestor("humanAxxxx","lemurAxxxx").get_leaves()) == 2
    assert len(T.get_common_ancestor("humanBxxxx","lemurBxxxx").get_leaves()) == 2
    assert len(T.get_common_ancestor("mouseAxxxx","ratAxxxxxx").get_leaves()) == 2
    assert len(T.get_common_ancestor("mouseBxxxx","ratBxxxxxx").get_leaves()) == 2
    assert len(T.get_common_ancestor("humanAxxxx","lemurAxxxx","mouseAxxxx","ratAxxxxxx").get_leaves()) == 4
    assert len(T.get_common_ancestor("humanBxxxx","lemurBxxxx","mouseBxxxx","ratBxxxxxx").get_leaves()) == 4
    assert len(T.get_common_ancestor("humanAxxxx","humanBxxxx").get_leaves()) == 8

    # Assert correct events inferred
    assert len(events) == 2
    assert events["D"] == 1
    assert events["S"] == 6
