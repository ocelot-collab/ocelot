from ocelot.cpbd.elements import Quadrupole


def test_quadrupole_k1l_getter():
    l = 2
    k1 = 0.3
    quad = Quadrupole(l=l, k1=k1)
    assert quad.k1l == l * k1

def test_quadrupole_k1l_setter():
    l = 2
    k1 = 0.3
    quad = Quadrupole(l=l, k1=k1)

    quad.k1l = 6
    assert quad.k1 == 3

def test_quadrupole_k2l_getter():
    l = 2
    k2 = 0.3
    quad = Quadrupole(l=l, k2=k2)
    assert quad.k2l == l * k2

def test_quadrupole_k2l_setter():
    l = 2.
    k2 = 0.3
    quad = Quadrupole(l=l, k2=k2)

    quad.k2l = 6
    assert quad.k2 == 3
    
