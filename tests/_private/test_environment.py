
import pytest
from unittest import mock

import os

from topiary._private.environment import load_env_variable
from topiary._private.check import check_float

def test_load_env_variable():

    with mock.patch.dict(os.environ, {'TEST1': '1.0', 'TEST2': 'one'}):

        # Just grab value without altering it
        value = load_env_variable("TEST1")
        assert value == "1.0"

        # Just grab value without altering it
        value = load_env_variable("TEST2")
        assert value == "one"

        # Just grab value without altering it
        value = load_env_variable("TEST3")
        assert value is None

        # Make sure check function parses value
        value = load_env_variable("TEST1",
                                  check_function=check_float)
        assert value == 1.0

        # Make sure check_function is being applied to check value
        with pytest.raises(ValueError):
            value = load_env_variable("TEST1",
                                      check_function=check_float,
                                      check_function_kwargs={"minimum_allowed":100.0})

        # Make sure check_function is being applied to check value
        with pytest.raises(ValueError):
            value = load_env_variable("TEST2",
                                      check_function=check_float)
        

                                    



