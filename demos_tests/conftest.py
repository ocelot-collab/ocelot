import os
import pytest

from params import *

def pytest_addoption(parser):
    parser.addoption('--update', action='store', default=None, help='used only for update reference values')


def pytest_collection_modifyitems(config, items):
    
    if config.getoption('--update') is not None:
        skip = pytest.mark.skip(reason='not used during update')
        for item in items:
            if 'update' not in item.keywords:
                item.add_marker(skip)
    else:
        skip = pytest.mark.skip(reason='used only for update')
        for item in items:
            if 'update' in item.keywords:
                item.add_marker(skip)


@pytest.fixture
def cmdopt(request):
    return request.config.getoption("--update")


# file with test results
pytest.TEST_RESULTS_FILE = os.path.dirname(os.path.abspath(__file__)) + '/test_results.dat'
f = open(pytest.TEST_RESULTS_FILE, 'w')
f.close()
