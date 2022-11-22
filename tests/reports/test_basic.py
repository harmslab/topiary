
import pytest

import topiary
from topiary.reports import session_to_html
from topiary.reports import df_to_table

import pandas as pd

def test_session_to_html():
    pass

def test_df_to_table():
    
    bad_input = [1,[],{},"test"]
    for b in bad_input:
        with pytest.raises(ValueError):
            df_to_table(b)

    df = pd.DataFrame({"test":[1,2,3],"this":["A","B","C"]})
    html = df_to_table(df)
    assert issubclass(type(html),str)
    assert html.startswith('<div class="table-responsive">')
    

    