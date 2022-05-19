
import pytest
import topiary


def test_create_pipeline_dict():

    pipe_dict = topiary._private.create_pipeline_dict()

    key_list = ["name","species","sequence","uid","keep"]

    assert set(pipe_dict.keys()) == set(key_list)

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
