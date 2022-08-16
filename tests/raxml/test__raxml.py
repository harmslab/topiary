import pytest

import topiary
from topiary.raxml._raxml import run_raxml
from topiary.raxml import RAXML_BINARY
from topiary._private import Supervisor

import os
import copy
import shutil

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_run_raxml(tiny_phylo,tmpdir):

    def _read_log_file(out_dir):

        f = open(os.path.join(out_dir,"test-log.sh"))
        out = f.read()
        f.close()

        return out.strip()

    df = tiny_phylo["final-output/dataframe.csv"]
    gene_tree = tiny_phylo["final-output/gene-tree.newick"]
    species_tree = tiny_phylo["final-output/species-tree.newick"]
    reconciled_tree = tiny_phylo["final-output/reconciled-tree.newick"]
    alignment = tiny_phylo["final-output/alignment.phy"]
    raxml_binary = shutil.which(RAXML_BINARY)

    current_dir = os.getcwd()
    os.chdir(tmpdir)

    # Real functionality is tested for specific phylogenetic situations via
    # things like test_raxml_tree, test_bootstrap, etc. Those functions rely
    # on run_raxml to correctly build raxml-ng commands, so focus on tests
    # to make sure args are correctly constructed. Rather than actually
    # running raxml, mostly write to test-log.sh because raxml will complain
    # about a lot of these command combos

    kwargs_template = {"algorithm":None,
                       "alignment_file":None,
                       "tree_file":None,
                       "model":None,
                       "seed":None,
                       "log_to_stdout":True,
                       "suppress_output":False,
                       "other_args":None,
                       "other_files":None,
                       "write_to_script":"test-log.sh",
                       "supervisor":None,
                       "num_threads":1,
                       "raxml_binary":"raxml-ng"}

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test0"
    kwargs["write_to_script"] = 1.0
    with pytest.raises(ValueError):
        run_raxml(**kwargs)

    # Make sure supervisor catches error
    sv = Supervisor()
    sv.create_calc_dir("test0.0","test_calc")
    assert sv.run_parameters["calc_status"] == "running"
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = os.path.join("test0.0","inside")
    kwargs["write_to_script"] = 1.0
    kwargs["supervisor"] = sv
    with pytest.raises(ValueError):
        run_raxml(**kwargs)
    assert sv.run_parameters["calc_status"] == "crashed"

    # Send in existing directory
    os.mkdir("test0.1")
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test0.1"
    with pytest.raises(FileExistsError):
        run_raxml(**kwargs)

    # Make sure supervisor catches error. (This throws error because
    # sv.create_calc_dir already created directory 0.2)
    sv = Supervisor()
    sv.create_calc_dir("test0.2","test_calc")
    assert sv.run_parameters["calc_status"] == "running"
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test0.2"
    kwargs["supervisor"] = sv
    with pytest.raises(FileExistsError):
        run_raxml(**kwargs)
    assert sv.run_parameters["calc_status"] == "crashed"

    # Send nothing in, get nothing out
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test1"
    run_raxml(**kwargs)
    output = _read_log_file("test1")
    assert output == f"{raxml_binary} --threads 1"

    # test binary. both tests for bad binary catch and that function is
    # actually reading binary argument
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test1.1"
    kwargs["raxml_binary"] = "raxml-ng-not-really-there"
    with pytest.raises(FileNotFoundError):
        run_raxml(**kwargs)

    # Make sure this error is caught by a supervisor if present
    sv = Supervisor()
    sv.create_calc_dir("test1.2","test_ml")
    assert sv.run_parameters["calc_status"] == "running"
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = os.path.join("test1.2","inside")
    kwargs["raxml_binary"] = "raxml-ng-not-really-there"
    kwargs["supervisor"] = sv
    with pytest.raises(FileNotFoundError):
        run_raxml(**kwargs)
    assert sv.run_parameters["calc_status"] == "crashed"

    # Send in alignment only
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test2"
    kwargs["alignment_file"] = alignment
    run_raxml(**kwargs)
    output = _read_log_file("test2")
    assert output == f"{raxml_binary} --msa alignment.phy --threads 1"
    assert os.path.isfile(os.path.join("test2","alignment.phy"))

    # Send in tree only
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test3"
    kwargs["tree_file"] = gene_tree
    run_raxml(**kwargs)
    output = _read_log_file("test3")
    assert output == f"{raxml_binary} --tree tree.newick --threads 1"
    assert os.path.isfile(os.path.join("test3","tree.newick"))

    # Send in other_files only
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test4"
    kwargs["other_files"] = [df,species_tree]
    run_raxml(**kwargs)
    output = _read_log_file("test4")
    assert output == f"{raxml_binary} --threads 1"
    assert os.path.isfile(os.path.join("test4","dataframe.csv"))
    assert os.path.isfile(os.path.join("test4","species-tree.newick"))

    # Send in algorithm only
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test5"
    kwargs["algorithm"] = "--somealgorithm"
    run_raxml(**kwargs)
    output = _read_log_file("test5")
    assert output == f"{raxml_binary} --somealgorithm --threads 1"

    # Send in model only
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test6"
    kwargs["model"] = "PRETEND_MODEL"
    run_raxml(**kwargs)
    output = _read_log_file("test6")
    assert output == f"{raxml_binary} --model PRETEND_MODEL --threads 1"

    # Send in algorithm, alignment, tree, model
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test7"
    kwargs["algorithm"] = "--somealgorithm"
    kwargs["alignment_file"] = alignment
    kwargs["tree_file"] = gene_tree
    kwargs["model"] = "PRETEND_MODEL"
    run_raxml(**kwargs)
    output = _read_log_file("test7")
    assert output == f"{raxml_binary} --somealgorithm --msa alignment.phy --tree tree.newick --model PRETEND_MODEL --threads 1"

    assert os.path.isfile(os.path.join("test7","alignment.phy"))
    assert os.path.isfile(os.path.join("test7","tree.newick"))

    # implicitly tested seed = None case above -- treated like False

    # Make sure a seed is being generated properly with True set
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test8"
    kwargs["seed"] = True
    run_raxml(**kwargs)
    output = _read_log_file("test8")
    assert output.startswith(f"{raxml_binary} --seed ")
    assert output[-11:] == "--threads 1"
    assert isinstance(int(output.split()[2]),int)

    # Make sure a seed is (not) being generated properly with False set
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test9"
    kwargs["seed"] = False
    run_raxml(**kwargs)
    output = _read_log_file("test9")
    assert output == f"{raxml_binary} --threads 1"

    # Make sure a seed is being generated properly with actual value set
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test10"
    kwargs["seed"] = 12345
    run_raxml(**kwargs)
    output = _read_log_file("test10")
    assert output == f"{raxml_binary} --seed 12345 --threads 1"

    # Make sure a seed is being generated set with integer-as-string sent in
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test11"
    kwargs["seed"] = "12345"
    run_raxml(**kwargs)
    output = _read_log_file("test11")
    assert output == f"{raxml_binary} --seed 12345 --threads 1"

    # Make sure throws error on stupid seed
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test12"
    kwargs["seed"] = "not_an_integer"
    with pytest.raises(ValueError):
        run_raxml(**kwargs)

    # Make sure supervisor catches error. (This throws error because
    # sv.create_calc_dir already created directory 0.2)
    sv = Supervisor()
    sv.create_calc_dir("test12.1","test_calc")
    assert sv.run_parameters["calc_status"] == "running"
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = os.path.join("test12.1","inside")
    kwargs["supervisor"] = sv
    kwargs["seed"] = "not_an_integer"
    with pytest.raises(ValueError):
        run_raxml(**kwargs)
    assert sv.run_parameters["calc_status"] == "crashed"

    # Auto generate threads
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test13"
    kwargs["num_threads"] = -1
    run_raxml(**kwargs)
    output = _read_log_file("test13")
    assert output.startswith(f"{raxml_binary} --threads ")
    assert int(output.split()[-1]) > 0

    # Specific number of threads
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test14"
    kwargs["num_threads"] = 2
    run_raxml(**kwargs)
    output = _read_log_file("test14")
    assert output == f"{raxml_binary} --threads 2"

    # Bad number of threads (this is tested by get_num_threads call, but
    # pretty important in this context; so check)
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test15"
    kwargs["num_threads"] = "not_a_number"
    with pytest.raises(ValueError):
        run_raxml(**kwargs)

    sv = Supervisor()
    sv.create_calc_dir("test15.1","test_calc")
    assert sv.run_parameters["calc_status"] == "running"
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = os.path.join("test15.1","inside")
    kwargs["supervisor"] = sv
    kwargs["num_threads"] = "not_a_number"
    with pytest.raises(ValueError):
        run_raxml(**kwargs)
    assert sv.run_parameters["calc_status"] == "crashed"

    # Test cross-talk between --all and --search flags and --num_threads
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test16"
    kwargs["algorithm"] = "--all"
    kwargs["num_threads"] = 2
    run_raxml(**kwargs)
    output = _read_log_file("test16")
    assert output == f"{raxml_binary} --all --threads auto" + "{2}"

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test17"
    kwargs["algorithm"] = "--search"
    kwargs["num_threads"] = 2
    run_raxml(**kwargs)
    output = _read_log_file("test17")
    assert output == f"{raxml_binary} --search --threads auto" + "{2}"

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test18"
    kwargs["algorithm"] = "--flagwithoutcrosstalk"
    kwargs["num_threads"] = 2
    run_raxml(**kwargs)
    output = _read_log_file("test18")
    assert output == f"{raxml_binary} --flagwithoutcrosstalk --threads 2"

    # Custom args
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test19"
    kwargs["other_args"] = ["other","args"]
    run_raxml(**kwargs)
    output = _read_log_file("test19")
    assert output == f"{raxml_binary} --threads 1 other args"

    # Send in supervisor. run_raxml should start and stop the job, which will
    # appear in the run_parameters dictionary
    sv = Supervisor()
    sv.create_calc_dir("test20","test_ml")
    assert len(sv.run_parameters["events"]) == 0

    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = os.path.join("test20","working","test")
    kwargs["supervisor"] = sv
    run_raxml(**kwargs)
    output = _read_log_file(os.path.join("test20","working","test"))
    assert output == f"{raxml_binary} --threads 1"
    assert os.path.isdir(sv.calc_dir)

    assert len(sv.run_parameters["events"]) == 1
    description = sv.run_parameters["events"][0]["description"]
    local_directory = sv.run_parameters["events"][0]["local_directory"]
    start_time = sv.run_parameters["events"][0]["time"]
    cmd = sv.run_parameters["events"][0]["cmd"]
    num_threads = sv.run_parameters["events"][0]["num_threads"]

    assert cmd[0] == f"{raxml_binary}"
    assert cmd[1] == "--threads"
    assert cmd[2] == "1"

    assert num_threads == 1

    # send in bad command with supervisor. make sure that function catches
    # and writes to supervisor before final RunTimeError.
    sv = Supervisor()
    sv.create_calc_dir("test21","test_ml")
    assert sv.run_parameters["calc_status"] == "running"

    # This calc will fail because we did not send anything into raxml except
    # the number of threads...
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = os.path.join("test21","working","test")
    kwargs["supervisor"] = sv
    kwargs["write_to_script"] = None
    with pytest.raises(RuntimeError):
        run_raxml(**kwargs)

    assert len(sv.run_parameters["events"]) == 1
    assert sv.run_parameters["calc_status"] == "crashed"
    # should have returned to starting directory after crash deep in working
    # dir
    assert os.getcwd() == tmpdir

    # Do one real calculation
    kwargs = copy.deepcopy(kwargs_template)
    kwargs["run_directory"] = "test22"
    kwargs["algorithm"] = "--search"
    kwargs["alignment_file"] = alignment
    kwargs["model"] = "LG"
    kwargs["write_to_script"] = None
    kwargs["seed"] = 3688946479

    output = run_raxml(**kwargs)
    assert output == f"{raxml_binary} --search --msa alignment.phy --model LG --seed 3688946479 --threads auto" + "{1}"
    assert os.path.isfile(os.path.join("test22","alignment.phy.raxml.bestTree"))

    f = open(gene_tree)
    expected = f.read().strip()
    f.close()

    f = open(os.path.join("test22","alignment.phy.raxml.bestTree"))
    observed = f.read().strip()
    f.close()

    # This is an aggressive test that requires raxml output be identical with
    # the same seed across platforms. Hopefully true. (Might fail on other boxes
    # because I made test file using experimental arm64 pll library).
    assert expected == observed

    os.chdir(current_dir)
