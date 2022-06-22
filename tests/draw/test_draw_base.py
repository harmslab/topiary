
import pytest

import numpy as np
import pandas as pd

import topiary

def test_ColorMap():

    expected_results = {(0,0,0):"#000000",
                        (255,0,0):"#ff0000",
                        (0,255,0):"#00ff00",
                        (255,255,0):"#ffff00",
                        (0,0,255):"#0000ff",
                        (255,0,255):"#ff00ff",
                        (0,255,255):"#00ffff",
                        (255,255,255):"#ffffff"}

    # Basic check of functionality. With defaults 0 -> 1 maps white -> black
    cm = topiary.draw._base.ColorMap()
    assert cm.hex(0) == expected_results[(255,255,255)]
    assert cm.hex(1) == expected_results[(0,0,0)]


    # Make sure input spans work
    good_span = [(-1,1),(0,1),(-1000,1000),(-1000,-500),[-1000,-500],
                 np.array((-1000,-500))]
    for g in good_span:
        print(f"trying span {g}")
        cm = topiary.draw._base.ColorMap(span=g)
        assert cm.hex(g[0]) == expected_results[(255,255,255)]
        assert cm.hex(g[1]) == expected_results[(0,0,0)]

    bad_span = [(1,1),(1,-1),None,{"test":1,"this":3},
                pd.DataFrame({"test":[1],"this":[3]}),"test"]
    for b in bad_span:
        print(f"trying span {b}")
        with pytest.raises(ValueError):
            cm = topiary.draw._base.ColorMap(span=b)

    # Make sure input channels work
    good_channel = [(0,1),(1,0),(0.5,1),(1,0.5),(0.5,0.5),(1,1),(0,0)]
    for g in good_channel:
        print(f"trying channel {g}")
        cm = topiary.draw._base.ColorMap(R=g)
        assert cm._R_bottom == g[0]
        assert cm._R_top == g[1]
        assert cm._R_intercept == g[0]
        assert cm._R_slope == g[1] - g[0]
        cm = topiary.draw._base.ColorMap(G=g)
        assert cm._G_bottom == g[0]
        assert cm._G_top == g[1]
        assert cm._G_intercept == g[0]
        assert cm._G_slope == g[1] - g[0]
        cm = topiary.draw._base.ColorMap(B=g)
        assert cm._B_bottom == g[0]
        assert cm._B_top == g[1]
        assert cm._B_intercept == g[0]
        assert cm._B_slope == g[1] - g[0]

    bad_channel = [(-1,1),(1,-1),(0,2),(1,2),None,{"test":1,"this":3},
                   pd.DataFrame({"test":[1],"this":[3]}),"test"]
    for b in bad_channel:
        print(f"trying channel {b}")
        with pytest.raises(ValueError):
            cm = topiary.draw._base.ColorMap(R=b)
        with pytest.raises(ValueError):
            cm = topiary.draw._base.ColorMap(G=b)
        with pytest.raises(ValueError):
            cm = topiary.draw._base.ColorMap(B=b)

    # Test individual channels, up and down
    cm = topiary.draw._base.ColorMap(R=(0,1),G=(0,0),B=(0,0))
    assert cm.hex(0) == expected_results[(0,0,0)]
    assert cm.hex(1) == expected_results[(255,0,0)]

    cm = topiary.draw._base.ColorMap(R=(1,0),G=(0,0),B=(0,0))
    assert cm.hex(0) == expected_results[(255,0,0)]
    assert cm.hex(1) == expected_results[(0,0,0)]

    cm = topiary.draw._base.ColorMap(R=(0,0),G=(0,1),B=(0,0))
    assert cm.hex(0) == expected_results[(0,0,0)]
    assert cm.hex(1) == expected_results[(0,255,0)]

    cm = topiary.draw._base.ColorMap(R=(0,0),G=(1,0),B=(0,0))
    assert cm.hex(0) == expected_results[(0,255,0)]
    assert cm.hex(1) == expected_results[(0,0,0)]

    cm = topiary.draw._base.ColorMap(R=(0,0),G=(0,0),B=(0,1))
    assert cm.hex(0) == expected_results[(0,0,0)]
    assert cm.hex(1) == expected_results[(0,0,255)]

    cm = topiary.draw._base.ColorMap(R=(0,0),G=(0,0),B=(0,1))
    assert cm.hex(0) == expected_results[(0,0,0)]
    assert cm.hex(1) == expected_results[(0,0,255)]

    cm = topiary.draw._base.ColorMap(R=(0,1),G=(0,1),B=(0,0))
    assert cm.hex(0) == expected_results[(0,0,0)]
    assert cm.hex(1) == expected_results[(255,255,0)]

    cm = topiary.draw._base.ColorMap(R=(0,1),G=(0,0),B=(0,1))
    assert cm.hex(0) == expected_results[(0,0,0)]
    assert cm.hex(1) == expected_results[(255,0,255)]

    cm = topiary.draw._base.ColorMap(R=(0,0),G=(0,1),B=(0,1))
    assert cm.hex(0) == expected_results[(0,0,0)]
    assert cm.hex(1) == expected_results[(0,255,255)]

    cm = topiary.draw._base.ColorMap(R=(0,1),G=(0,1),B=(0,1))
    assert cm.hex(0) == expected_results[(0,0,0)]
    assert cm.hex(1) == expected_results[(255,255,255)]
