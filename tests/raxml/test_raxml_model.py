
import pytest

import topiary

from topiary.raxml.model import _generate_parsimony_tree
from topiary.raxml.model import _model_thread_function
from topiary.raxml.model import _parse_raxml_info_for_aic
from topiary.raxml.model import find_best_model

import pandas as pd

import os, json

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test__generate_parsimony_tree():

    pass

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test__model_thread_function():

    pass

@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test__parse_raxml_info_for_aic():

    pass


@pytest.mark.skipif(os.name == "nt",reason="cannot run on windows")
def test_find_best_model(tiny_phylo,tmpdir):

    df = tiny_phylo["initial-input/dataframe.csv"]
    current_dir = os.getcwd()
    os.chdir(tmpdir)

    model_matrices = ["LG","JTT"]
    model_rates = ["","G8"]
    model_freqs = ["","FO"]
    model_invariant = ["","IO"]
    calc_dir = "test_out"

    find_best_model(df,
                    model_matrices=model_matrices,
                    model_rates=model_rates,
                    model_freqs=model_freqs,
                    model_invariant=model_invariant,
                    calc_dir=calc_dir)

    # Make sure we're building models properly
    all_models = []
    for a in model_matrices:
        for b in model_rates:
            for c in model_freqs:
                for d in model_invariant:
                    model = [m for m in [a,b,c,d] if m != ""]
                    all_models.append("+".join(model))

    # Read output dataframe
    out_df = pd.read_csv(os.path.join("test_out","output","model-comparison.csv"))

    # Make sure it is the right length, sorted correctly, and has all models
    assert len(out_df) == len(all_models)
    assert out_df["p"].iloc[0] > out_df["p"].iloc[-1]
    assert len(set(all_models).intersection(set(out_df["model"]))) == len(all_models)

    # Make sure we can read sane json
    f = open(os.path.join("test_out","run_parameters.json"))
    out_json = json.load(f)
    f.close()

    # Make sure best model is recorded properly in output json
    best_model = out_df["model"].iloc[0]
    assert out_json["model"] == best_model
    assert out_json["calc_status"] == "complete"

    os.chdir(current_dir)
