
import pytest
from unittest import mock

from topiary._private.mpi import get_hosts
from topiary._private.mpi import get_num_slots
from topiary._private.mpi import check_mpi_configuration

from topiary.generax import GENERAX_BINARY

import os

@pytest.mark.run_generax
def test_get_hosts():

    hosts = get_hosts(1)
    assert len(hosts) == 1
    assert isinstance(hosts[0],str)

    hosts = get_hosts(2)
    assert len(hosts) == 2
    assert isinstance(hosts[0],str)
    assert isinstance(hosts[1],str)


@pytest.mark.run_generax
def test_get_num_slots():

    # This should work for any system with more tha one core. 
    num_slots = get_num_slots()
    assert num_slots > 1

    # Enforce a single slot with an environment variable
    with mock.patch.dict(os.environ, {'TOPIARY_MAX_SLOTS': '1'}):
        num_slots = get_num_slots()
        assert num_slots == 1


@pytest.mark.run_generax
def test_check_mpi_configuration():

    # This should run fine
    check_mpi_configuration(1)

    # This should run fine (assuming test environment has more than one slot)
    check_mpi_configuration(2)

    # An unlikely number of slots. If this test ever passes, we've reached the
    # singularity and this code does not matter anyway
    with pytest.raises(RuntimeError):
        check_mpi_configuration(1000000)
