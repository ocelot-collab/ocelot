import numpy as np
import pytest

from unit_tests.params import check_dict, check_matrix, check_value


@pytest.mark.parametrize(
    'value, reference',
    [(4.2603605111289117e-10, 4.260361033148104e-10),
     (-1.4462257100222864e-12, -1.4462444391377672e-12)],
)
def test_absolute_allowance_accepts_roundoff_near_zero(value, reference):
    assert check_value(value, reference, tolerance=1e-7) is not None
    assert check_value(value, reference, tolerance=1e-7,
                       absolute_tolerance=1e-15) is None


@pytest.mark.parametrize('value, reference', [(2e-12, 1e-12), (1.001, 1.0)])
def test_absolute_allowance_still_rejects_real_changes(value, reference):
    assert check_value(value, reference, tolerance=1e-7,
                       absolute_tolerance=1e-15) is not None


def test_absolute_allowance_propagates_through_nested_reference_data():
    actual = [{'scalar': 1e-12 + 2e-17,
               'array': [1e-12 + 2e-17],
               'nested': {'value': 1e-12 + 2e-17}}]
    reference = [{'scalar': 1e-12,
                  'array': [1e-12],
                  'nested': {'value': 1e-12}}]
    assert all(r is not None for r in check_dict(actual, reference, tolerance=1e-7))
    assert check_dict(actual, reference, tolerance=1e-7,
                      absolute_tolerance=1e-15) == [None, None, None]
    assert check_matrix(np.array([actual[0]['scalar']]), np.array([1e-12]),
                        tolerance=1e-7, absolute_tolerance=1e-15) == [None]


def test_absolute_mode_keeps_its_existing_threshold():
    assert check_value(2e-12, 1e-12, tolerance=1e-13, tolerance_type='absolute',
                       absolute_tolerance=1.0) is not None
