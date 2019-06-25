"""Test parameters description file"""

import pytest
import numpy as np

from ocelot import *
from ocelot.cpbd.beam import generate_parray





@pytest.fixture(scope='function')
def p_array():

    np.random.seed(11)

    # generate beam file
    sigma_x = 0.000121407185261
    sigma_px = 1.80989470506e-05
    sigma_y = 0.000165584800564
    sigma_py = 4.00994225888e-05

    p_array = generate_parray(sigma_x=sigma_x, sigma_px=sigma_px, sigma_y=sigma_y, sigma_py=sigma_py,
                              sigma_tau=1.30190131e-04, sigma_p=3.09815718e-04, chirp=0.002, charge=0.5e-9,
                              nparticles=100, energy=0.13)


    return p_array