import numpy as np

from ocelot.cpbd.elements import Bend, Cavity, Drift, Matrix, Quadrupole
from ocelot.cpbd.latticeIO import LatticeIO
from ocelot.cpbd.magnetic_lattice import MagneticLattice


def _load_namespace(lines):
    namespace = {}
    exec("".join(lines), namespace)
    return namespace


def test_get_elements_returns_unique_sequence_elements_in_order():
    drift = Drift(l=0.2, eid="D.1")
    quad = Quadrupole(l=0.3, k1=1.4, eid="Q.1")
    lattice = MagneticLattice((drift, quad, drift))

    elements = LatticeIO._get_elements(lattice)

    assert elements == [drift, quad]
    assert drift.name == "d_1"
    assert quad.name == "q_1"


def test_lat2input_round_trips_bend_and_cavity_state():
    bend = Bend(l=1.2, angle=0.03, k1=-0.2, e1=0.01, e2=-0.02, eid="B1")
    cavity = Cavity(
        l=0.5,
        v=0.02,
        phi=30.0,
        freq=1.3e9,
        vx_up=-1.1e-5 + 2.2e-6j,
        vy_up=3.3e-5 - 4.4e-6j,
        vxx_up=5.5e-4 + 6.6e-5j,
        vxy_up=-7.7e-4 + 8.8e-5j,
        vx_down=9.9e-6 - 1.1e-6j,
        vy_down=-2.2e-5 + 3.3e-6j,
        vxx_down=4.4e-4 - 5.5e-5j,
        vxy_down=-6.6e-4 - 7.7e-5j,
        eid="C1",
    )
    lattice = MagneticLattice((bend, cavity))

    namespace = _load_namespace(LatticeIO.lat2input(lattice))
    rebuilt = MagneticLattice(namespace["cell"], method=lattice.method)
    rebuilt_bend, rebuilt_cavity = rebuilt.sequence

    assert isinstance(rebuilt_bend, Bend)
    assert rebuilt_bend.id == bend.id
    assert rebuilt_bend.l == bend.l
    assert rebuilt_bend.angle == bend.angle
    assert rebuilt_bend.k1 == bend.k1
    assert rebuilt_bend.e1 == bend.e1
    assert rebuilt_bend.e2 == bend.e2

    assert isinstance(rebuilt_cavity, Cavity)
    assert rebuilt_cavity.id == cavity.id
    assert rebuilt_cavity.l == cavity.l
    assert rebuilt_cavity.v == cavity.v
    assert rebuilt_cavity.phi == cavity.phi
    assert rebuilt_cavity.freq == cavity.freq
    assert rebuilt_cavity.vx_up == cavity.vx_up
    assert rebuilt_cavity.vy_up == cavity.vy_up
    assert rebuilt_cavity.vxx_up == cavity.vxx_up
    assert rebuilt_cavity.vxy_up == cavity.vxy_up
    assert rebuilt_cavity.vx_down == cavity.vx_down
    assert rebuilt_cavity.vy_down == cavity.vy_down
    assert rebuilt_cavity.vxx_down == cavity.vxx_down
    assert rebuilt_cavity.vxy_down == cavity.vxy_down


def test_lat2input_round_trips_matrix_arrays():
    matrix = Matrix(l=0.4, delta_e=0.01, eid="M1")
    matrix.r[0, 0] = 1.1
    matrix.r[1, 5] = -0.02
    matrix.t[0, 0, 0] = 0.3
    matrix.t[4, 5, 5] = -0.4
    matrix.b[0, 0] = 0.05
    matrix.b[3, 0] = -0.07
    lattice = MagneticLattice((matrix,))

    namespace = _load_namespace(LatticeIO.lat2input(lattice))
    rebuilt_matrix = MagneticLattice(namespace["cell"], method=lattice.method).sequence[0]

    assert isinstance(rebuilt_matrix, Matrix)
    assert rebuilt_matrix.id == matrix.id
    assert rebuilt_matrix.l == matrix.l
    assert rebuilt_matrix.delta_e == matrix.delta_e
    np.testing.assert_allclose(rebuilt_matrix.r, matrix.r)
    np.testing.assert_allclose(rebuilt_matrix.t, matrix.t)
    np.testing.assert_allclose(rebuilt_matrix.b, matrix.b)
