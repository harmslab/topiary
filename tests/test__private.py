
import pytest
import topiary
from topiary._private import wrap_function as wf

import numpy as np

import os

def test_generate_uid():

    uid = topiary._private.generate_uid()
    assert type(uid) is str
    assert len(uid) == 10

    uid = topiary._private.generate_uid(2)
    assert type(uid) is list
    assert len(uid) == 2
    assert len(list(set(uid))) == 2
    for i in range(2):
        assert type(uid[i]) is str
        assert len(uid[i]) == 10


def test_wrap_function(tmpdir):

    def test_fcn(arg1,
                 arg2=None,
                 arg3=["test","this"],
                 arg4=True,
                 arg5=False,
                 arg6=1.0,
                 arg7="test",
                 arg8=[1.0,2.0]):

        return arg1, arg2, arg3, arg4, arg5, arg6, arg7

    with pytest.raises(TypeError):
        out = wf()

    out = wf(test_fcn,argv=["stupid"])
    assert out.__dict__["arg2"] is None
    assert np.array_equal(out.__dict__["arg3"],["test","this"])
    assert out.__dict__["arg4"] is True
    assert out.__dict__["arg5"] is False
    assert np.isclose(out.__dict__["arg6"],1.0)
    assert out.__dict__["arg7"] == "test"
    assert np.array_equal(out.__dict__["arg8"],[1,2])

    # Don't send in useful first arg to argv. Should throw error because arg
    # is required.
    with pytest.raises(SystemExit):
        out = wf(test_fcn,argv=["--arg2","numb"])

    # Check basic parsing of the args. Pass in arg2 and make sure arg3-7 are
    # unchanged.
    out = wf(test_fcn,argv=["stupid","--arg2","numb"])
    assert out.__dict__["arg2"] == "numb"
    assert np.array_equal(out.__dict__["arg3"],["test","this"])
    assert out.__dict__["arg4"] is True
    assert out.__dict__["arg5"] is False
    assert np.isclose(out.__dict__["arg6"],1.0)
    assert out.__dict__["arg7"] == "test"
    assert np.array_equal(out.__dict__["arg8"],[1,2])

    # Make sure arg3-7 are parsed properly when sent in
    out = wf(test_fcn,argv=["stupid","--arg3","bob"])
    #assert np.array_equal(out.__dict__["arg3"],list("bob"))

    out = wf(test_fcn,argv=["stupid","--arg4"])
    assert out.__dict__["arg4"] == False

    out = wf(test_fcn,argv=["stupid","--arg5"])
    assert out.__dict__["arg5"] == True

    out = wf(test_fcn,argv=["stupid","--arg6","5"])
    assert out.__dict__["arg6"] == 5.0

    out = wf(test_fcn,argv=["stupid","--arg6","-5"])
    assert out.__dict__["arg6"] == -5.0

    out = wf(test_fcn,argv=["stupid","--arg7","string"])
    assert out.__dict__["arg7"] == "string"

    # Make sure arg type checking working (really only matters for float)
    with pytest.raises(SystemExit):
        out = wf(test_fcn,argv=["stupid","--arg6","string"])

    optional_arg_types = {"arg2":str}
    out = wf(test_fcn,argv=["stupid","--arg2","numb"],optional_arg_types=optional_arg_types)
    assert out.__dict__["arg2"] == "numb"

    optional_arg_types = {"arg2":float}
    with pytest.raises(SystemExit):
        out = wf(test_fcn,argv=["stupid","--arg2","numb"],optional_arg_types=optional_arg_types)

    optional_arg_types = {"arg2":float}
    out = wf(test_fcn,argv=["stupid","--arg2","1.0"],optional_arg_types=optional_arg_types)
    assert out.__dict__["arg2"] == 1.0

    # Make sure list parsing works
    # Make sure list parsing works
    out = wf(test_fcn,argv=["stupid","--arg8","3"])
    assert np.array_equal(out.__dict__["arg8"],[3])

    out = wf(test_fcn,argv=["stupid","--arg8","3","4"])
    assert np.array_equal(out.__dict__["arg8"],[3,4])

    out = wf(test_fcn,argv=["stupid","--arg8","3","4","5"])
    assert np.array_equal(out.__dict__["arg8"],[3,4,5])

    # Test file read for value
    test_file = os.path.join(tmpdir,'test.txt')
    f = open(test_file,'w')
    f.write("3\n4\n5\n")
    f.close()

    out = wf(test_fcn,argv=["stupid","--arg8",test_file])
    assert np.array_equal(out.__dict__["arg8"],[3,4,5])

    # Test comment/blank line parsing
    test_file = os.path.join(tmpdir,'test.txt')
    f = open(test_file,'w')
    f.write("#3\n4\n5\n\n\n")
    f.close()

    out = wf(test_fcn,argv=["stupid","--arg8",test_file])
    assert np.array_equal(out.__dict__["arg8"],[4,5])

    # Test bad file  (with non-float values)
    test_file = os.path.join(tmpdir,'test.txt')
    f = open(test_file,'w')
    f.write("A\nB\n5\n\n\n")
    f.close()
    with pytest.raises(ValueError):
        out = wf(test_fcn,argv=["stupid","--arg8",test_file])
