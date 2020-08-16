"""Test parameters description file"""

import pytest

from ocelot import *


"""pytest fixtures definition"""
    
@pytest.fixture(scope='module')
def method():

    m = MethodTM()
    m.global_method = SecondTM
    
    return m

    
@pytest.fixture(scope='function')
def tws0():
    
    tws0 = Twiss()
    tws0.beta_x = 13.85
    tws0.beta_y = 13.85
    tws0.alpha_x = -1.695
    tws0.alpha_y = -1.692
    tws0.E = 0.155
    
    return tws0
