
import pytest
import topiary
from topiary._private import wrap_function as wf

import numpy as np

import os


def test_wrap_function(tmpdir):

    def test_fcn(arg1,
                 arg2=None,
                 arg3=["test","this"],
                 arg4=True,
                 arg5=False,
                 arg6=1.0,
                 arg7="test",
                 arg8=[1.0,2.0]):
        """
        test function.
        """

        return arg1, arg2, arg3, arg4, arg5, arg6, arg7

    with pytest.raises(TypeError):
        ret, out = wf()

    ret, out = wf(test_fcn,argv=["stupid"])
    assert out.__dict__["arg1"] == "stupid"
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
        ret, out = wf(test_fcn,argv=["--arg2","numb"])

    # Check basic parsing of the args. Pass in arg2 and make sure arg3-7 are
    # unchanged.
    ret, out = wf(test_fcn,argv=["stupid","--arg2","numb"])
    assert out.__dict__["arg1"] == "stupid"
    assert out.__dict__["arg2"] == "numb"
    assert np.array_equal(out.__dict__["arg3"],["test","this"])
    assert out.__dict__["arg4"] is True
    assert out.__dict__["arg5"] is False
    assert np.isclose(out.__dict__["arg6"],1.0)
    assert out.__dict__["arg7"] == "test"
    assert np.array_equal(out.__dict__["arg8"],[1,2])

    # Make sure arg3-7 are parsed properly when sent in
    ret, out = wf(test_fcn,argv=["stupid","--arg3","bob"])
    #assert np.array_equal(out.__dict__["arg3"],list("bob"))

    ret, out = wf(test_fcn,argv=["stupid","--arg4"])
    assert out.__dict__["arg4"] == False

    ret, out = wf(test_fcn,argv=["stupid","--arg5"])
    assert out.__dict__["arg5"] == True

    ret, out = wf(test_fcn,argv=["stupid","--arg6","5"])
    assert out.__dict__["arg6"] == 5.0

    ret, out = wf(test_fcn,argv=["stupid","--arg6","-5"])
    assert out.__dict__["arg6"] == -5.0

    ret, out = wf(test_fcn,argv=["stupid","--arg7","string"])
    assert out.__dict__["arg7"] == "string"

    # Make sure arg type checking working (really only matters for float)
    with pytest.raises(SystemExit):
        ret, out = wf(test_fcn,argv=["stupid","--arg6","string"])

    optional_arg_types = {"arg2":str}
    ret, out = wf(test_fcn,argv=["stupid","--arg2","numb"],optional_arg_types=optional_arg_types)
    assert out.__dict__["arg2"] == "numb"

    optional_arg_types = {"arg2":float}
    with pytest.raises(SystemExit):
        ret, out = wf(test_fcn,argv=["stupid","--arg2","numb"],optional_arg_types=optional_arg_types)

    optional_arg_types = {"arg2":float}
    ret, out = wf(test_fcn,argv=["stupid","--arg2","1.0"],optional_arg_types=optional_arg_types)
    assert out.__dict__["arg2"] == 1.0

    # Make sure list parsing works
    ret, out = wf(test_fcn,argv=["stupid","--arg8","3"])
    assert np.array_equal(out.__dict__["arg8"],[3])

    ret, out = wf(test_fcn,argv=["stupid","--arg8","3","4"])
    assert np.array_equal(out.__dict__["arg8"],[3,4])

    ret, out = wf(test_fcn,argv=["stupid","--arg8","3","4","5"])
    assert np.array_equal(out.__dict__["arg8"],[3,4,5])

    # Test file read for value
    test_file = os.path.join(tmpdir,'test.txt')
    f = open(test_file,'w')
    f.write("3\n4\n5\n")
    f.close()

    ret, out = wf(test_fcn,argv=["stupid","--arg8",test_file])
    assert np.array_equal(out.__dict__["arg8"],[3,4,5])

    # Test comment/blank line parsing
    test_file = os.path.join(tmpdir,'test.txt')
    f = open(test_file,'w')
    f.write("#3\n4\n5\n\n\n")
    f.close()

    ret, out = wf(test_fcn,argv=["stupid","--arg8",test_file])
    assert np.array_equal(out.__dict__["arg8"],[4,5])

    # Test bad file  (with non-float values)
    test_file = os.path.join(tmpdir,'test.txt')
    f = open(test_file,'w')
    f.write("A\nB\n5\n\n\n")
    f.close()
    with pytest.raises(ValueError):
        ret, out = wf(test_fcn,argv=["stupid","--arg8",test_file])

    # Extra arguments

    # Do not pass in required extra argument
    with pytest.raises(SystemExit):
        ret, out = wf(test_fcn,argv=["stupid"],extra_args=[("extra",{"type":str})])

    # string required
    ret, out = wf(test_fcn,argv=["stupid","rocket"],
                  extra_args=[("extra",{"type":str})])
    assert out.__dict__["extra"] == "rocket"

    # int required
    ret, out = wf(test_fcn,argv=["stupid","5"],
                  extra_args=[("extra",{"type":int})])
    assert out.__dict__["extra"] == 5


    # string optional
    ret, out = wf(test_fcn,argv=["stupid"],
                  extra_args=[("--extra",{"type":str,"default":"xnay"})])
    assert out.__dict__["extra"] == "xnay"

    ret, out = wf(test_fcn,argv=["stupid","--extra","rocket"],
                  extra_args=[("--extra",{"type":str})])
    assert out.__dict__["extra"] == "rocket"

    # int optional
    ret, out = wf(test_fcn,argv=["stupid","--extra","5"],
                  extra_args=[("--extra",{"type":int})])
    assert out.__dict__["extra"] == 5

    ret, out = wf(test_fcn,argv=["stupid"],
                  extra_args=[("--extra",{"type":int,"default":-5})])
    assert out.__dict__["extra"] == -5

    # bool optional
    ret, out = wf(test_fcn,argv=["stupid","--extra"],
                  extra_args=[("--extra",{"action":"store_false"})])
    assert out.__dict__["extra"] == False

    ret, out = wf(test_fcn,argv=["stupid"],
                  extra_args=[("--extra",{"action":"store_false"})])
    assert out.__dict__["extra"] == True

    ret, out = wf(test_fcn,argv=["stupid","--extra"],
                  extra_args=[("--extra",{"action":"store_true"})])
    assert out.__dict__["extra"] == True

    ret, out = wf(test_fcn,argv=["stupid"],
                  extra_args=[("--extra",{"action":"store_true"})])
    assert out.__dict__["extra"] == False
