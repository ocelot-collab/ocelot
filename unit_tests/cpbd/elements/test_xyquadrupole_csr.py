"""CSR geometry follows the selected XYQuadrupole reference model."""
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from ocelot.cpbd.csr import CSR
from ocelot.cpbd.elements import Bend, Drift, XYQuadrupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice


def trajectory(elements, rk=False):
    lattice = MagneticLattice(elements)
    csr = CSR()
    csr.indx0, csr.indx1 = 0, len(elements)-1
    csr.energy = 0.001
    csr.rk_traj = rk
    csr.traj_step = 0.005
    return lattice, csr.prepare(lattice)


@pytest.mark.parametrize('k1', [-1., 0., 1.])
@pytest.mark.parametrize('tilt', [0., .4])
def test_csr_analytical_xyquad_matches_survey_with_preceding_bend(k1, tilt):
    q = XYQuadrupole(l=.8, k1=k1, x_offs=.01, y_offs=.03, tilt=tilt)
    lattice, traj = trajectory([Bend(l=.2, angle=.03, tilt=.3), q, Drift(l=.1)])
    _, survey = lattice.survey()
    np.testing.assert_allclose(traj[1:4, -1], survey[-1]['r_end'], atol=2e-12, rtol=0)
    np.testing.assert_allclose(traj[4:7, -1], survey[-1]['W'][:, 2], atol=2e-12, rtol=0)
    np.testing.assert_allclose(np.linalg.norm(traj[4:7], axis=0), 1, atol=2e-12)
    assert traj[0, -1] == pytest.approx(lattice.totalLen)
    assert np.all(np.diff(traj[0]) > 0)


@pytest.mark.parametrize('tilt', [0., .4])
def test_csr_rk_xyquad_matches_independent_lorentz_reference(tilt):
    q = XYQuadrupole(l=.8, k1=1.2, x_offs=.01, y_offs=.03, tilt=tilt)
    _, traj = trajectory([q], rk=True)
    c, s = np.cos(tilt), np.sin(tilt)
    rotation = np.array([[c,-s,0], [s,c,0], [0,0,1]])

    def rhs(arc, state):
        x, y, z = rotation.T @ state[:3]
        field = rotation @ np.array([q.k1*(y+q.y_offs), q.k1*(x+q.x_offs), 0])
        return np.r_[state[3:], np.cross(state[3:], field)]

    expected = solve_ivp(rhs, [0, q.l], [0,0,0,0,0,1], t_eval=traj[0],
                         method='DOP853', rtol=2e-13, atol=2e-15).y
    np.testing.assert_allclose(traj[1:], expected, atol=3e-12, rtol=0)
    # The curvature grows along the defocusing plane: it is not a uniform dipole.
    changes = np.linalg.norm(np.diff(traj[4:7], axis=1), axis=0)
    assert changes[-1] > changes[0]*1.2
