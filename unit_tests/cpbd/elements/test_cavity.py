from random import random
from ocelot.cpbd.elements import Cavity



def test_cavity_remove_coupler_kick():
    ckicks = dict(vx_up=random(), vy_up=random(), vxx_up=random(),
                  vxy_up=random(), vx_down=random(), vy_down=random(),
                  vxx_down=random(), vxy_down=random())

    cav = Cavity(l=3.0, phi=4, freq=3e6, **ckicks)

    cav.remove_coupler_kick()
    assert cav.vx_up == 0
    assert cav.vy_up == 0
    assert cav.vxx_up == 0
    assert cav.vxy_up == 0
    assert cav.vx_down == 0
    assert cav.vy_down == 0
    assert cav.vxx_down == 0
    assert cav.vxy_down == 0
