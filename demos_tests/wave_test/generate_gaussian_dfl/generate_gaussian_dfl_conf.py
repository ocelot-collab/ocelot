"""Test parameters description file"""

import pytest


@pytest.fixture(scope='module')
def args_array_ref():
    args_array_ref = [ret_args(1e-10, shape=(30, 30, 5), dgrid=(1e-3, 1e-3, None),
                               power_rms=(0.2e-3, 0.2e-3, 0.1e-6), power_center=(0, 0, None),
                               power_angle=(0, 0), power_waistpos=(0.1e-3, 0), wavelength=None,
                               zsep=1, freq_chirp=0, en_pulse=None, power=10e6, debug=1),
                      ret_args(1e-10, shape=(30, 30, 5), dgrid=(1e-3, 1e-3, None),
                               power_rms=(0.2e-3, 0.2e-3, 0.1e-6), power_center=(0, 0.1e-3, None),
                               power_angle=(0, 0), power_waistpos=(0, 0), wavelength=None,
                               zsep=1, freq_chirp=0, en_pulse=None, power=10e6, debug=1),
                      ret_args(1e-10, shape=(30, 30, 5), dgrid=(1e-3, 1e-3, None),
                               power_rms=(0.4e-3, 0.2e-3, 0.1e-6), power_center=(0, -0.1e-3, None),
                               power_angle=(0, 0), power_waistpos=(0, 0), wavelength=None,
                               zsep=1, freq_chirp=0, en_pulse=None, power=10e6, debug=1),
                      ret_args(1e-10, shape=(30, 30, 5), dgrid=(1e-3, 1e-3, None),
                               power_rms=(0.4e-3, 0.4e-3, 0.1e-6), power_center=(-0.2e-3, 0.5e-3, None),
                               power_angle=(0, 0), power_waistpos=(0, 0), wavelength=None,
                               zsep=1, freq_chirp=0, en_pulse=None, power=10e6, debug=1)]
    return args_array_ref


@pytest.fixture
def cmdopt(request):
    return request.config.getoption("--update")


def ret_args(*args, **kwargs):
    return (args, kwargs)
