"""Test parameters description file"""

import pytest


@pytest.fixture(scope='function')
def args_array_ref():
    args_array_ref = [ret_args(1e-10, shape=(21, 21, 21), dgrid=(1e-3, 1e-3, None),
                                power_rms=(0.1e-3, 0.1e-3, 0.01e-6), power_center=(0, 0, None),
                                power_angle=(0, 0), power_waistpos=(0.1e-3, 0), wavelength=None,
                                zsep=20, freq_chirp=0, en_pulse=None, power=10e6),
                      ret_args(1e-10, shape=(21, 21, 21), dgrid=(1e-3, 1e-3, 0.04e-6),
                                power_rms=(0.1e-3, 0.1e-3, 0.01e-6), power_center=(-0.25e-3, -0.25e-3, None),
                                power_angle=(0.5e-6, 0.5e-6), power_waistpos=(0.1e-3, 0), wavelength=None,
                                zsep=None, freq_chirp=0, en_pulse=None, power=10e6),
                      ret_args(1e-10, shape=(11, 11, 500), dgrid=(1e-3, 1e-3, None),
                                power_rms=(0.1e-3, 0.1e-3, 0.01e-6), power_center=(0, 0, None),
                                power_angle=(0, 0), power_waistpos=(0.1e-3, 0), wavelength=1.2e-10, 
                                zsep=1, freq_chirp=8e-4, en_pulse=1e-6, power=None),
                      # ret_args(1e-10, shape=(11, 11, 5), dgrid=(500e-6, 500e-6, None),
                                # power_rms=(100e-6, 100e-6, 2e-9), power_center=(0, 0, None),
                                # power_angle=(0, 0), power_waistpos=(0, 0), wavelength=None,
                                # zsep=20, freq_chirp=0, en_pulse=None, power=10e6),
                      # ret_args(1e-10, shape=(11, 11, 5), dgrid=(500e-6, 500e-6, 10e-9),
                                # power_rms=(100e-6, 100e-6, 2e-9), power_center=(-150e-6, -150e-6, None),
                                # power_angle=(0.5e-6, 0.5e-6), power_waistpos=(0, 0), wavelength=None,
                                # zsep=None, freq_chirp=0, en_pulse=None, power=10e6,
                      # ret_args(1e-10, shape=(3, 3, 20), dgrid=(500e-6, 500e-6, None),
                                # power_rms=(100e-6, 100e-6, 0.1e-9), power_center=(0, 0, None),
                                # power_angle=(0, 0), power_waistpos=(0, 0), wavelength=None,
                                # zsep=1, freq_chirp=0.1, en_pulse=None, power=10e6,
                               ]
    return args_array_ref


@pytest.fixture
def cmdopt(request):
    return request.config.getoption("--update")


def ret_args(*args, **kwargs):
    return (args, kwargs)
