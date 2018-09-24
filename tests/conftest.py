"""pytest configuration for kdotp-symmetry tests."""
# pylint: disable=unused-argument,redefined-outer-name,protected-access

import json
import pytest


@pytest.fixture
def test_name(request):
    """Returns module_name.function_name for a given test"""
    return request.module.__name__ + '/' + request._parent_request._pyfuncitem.name


@pytest.fixture
def compare_data(request, test_name, scope="session"):
    """Returns a function which either saves some data to a file or (if that file exists already) compares it to pre-existing data using a given comparison function."""

    def inner(compare_fct, data, tag=None):
        full_name = test_name + (tag or '')

        # get rid of json-specific quirks
        # store as string because I cannot add the decoder to the pytest cache
        data_str = json.dumps(data)
        data = json.loads(data_str)
        val = json.loads(request.config.cache.get(full_name, 'null'))

        if val is None:
            request.config.cache.set(full_name, data_str)
            raise ValueError('Reference data does not exist.')
        else:
            assert compare_fct(val, data)

    return inner


@pytest.fixture
def compare_equal(compare_data):
    """Returns a function which compares data for equality with pre-existing data."""
    return lambda data, tag=None: compare_data(lambda x, y: x == y, data, tag)
