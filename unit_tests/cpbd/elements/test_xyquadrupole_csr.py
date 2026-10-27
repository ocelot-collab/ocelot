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
    q = XYQuadrupole(l=.8, k1=k1, x_offs=.01, y_offs=.03, tilt=tilt, x_angle=.02, y_angle=-.01)
    lattice, traj = trajectory([Bend(l=.2, angle=.03, tilt=.3), q, Drift(l=.1)])
    _, survey = lattice.survey()
    np.testing.assert_allclose(traj[1:4, -1], survey[-1]['r_end'], atol=2e-12, rtol=0)
    np.testing.assert_allclose(traj[4:7, -1], survey[-1]['W'][:, 2], atol=2e-12, rtol=0)
    np.testing.assert_allclose(np.linalg.norm(traj[4:7], axis=0), 1, atol=2e-12)
    assert traj[0, -1] == pytest.approx(lattice.totalLen)
    assert np.all(np.diff(traj[0]) > 0)


@pytest.mark.parametrize('tilt', [0., .4])
@pytest.mark.parametrize('angles', [(0., 0.), (.02, -.01)])
def test_csr_rk_xyquad_matches_independent_lorentz_reference(tilt, angles):
    from scipy.spatial.transform import Rotation

    q = XYQuadrupole(l=.8, k1=1.2, x_offs=.01, y_offs=.03, tilt=tilt,
                     x_angle=angles[0], y_angle=angles[1])
    _, traj = trajectory([q], rk=True)
    c, s = np.cos(tilt), np.sin(tilt)
    rotation = np.array([[c,-s,0], [s,c,0], [0,0,1]])
    transverse = np.hypot(*angles)
    entrance = np.eye(3)
    if transverse:
        axis = np.array([-angles[1], angles[0], 0.]) / transverse
        entrance = Rotation.from_rotvec(np.arcsin(transverse)*axis).as_matrix()

    def rhs(arc, state):
        x, y, z = state[:3]
        field = np.array([q.k1*y, q.k1*x, 0])
        return np.r_[state[3:], np.cross(state[3:], field)]

    initial_position = np.array([q.x_offs, q.y_offs, 0.])
    solution = solve_ivp(rhs, [0, q.l], np.r_[initial_position, entrance[:, 2]], t_eval=traj[0],
                         method='DOP853', rtol=2e-13, atol=2e-15)
    assert solution.success
    expected = np.vstack((rotation @ entrance.T @ (solution.y[:3] - initial_position[:, None]),
                          rotation @ entrance.T @ solution.y[3:]))
    np.testing.assert_allclose(traj[1:], expected, atol=3e-12, rtol=0)
    # The curvature grows along the defocusing plane: it is not a uniform dipole.
    changes = np.linalg.norm(np.diff(traj[4:7], axis=1), axis=0)
    assert changes[-1] > changes[0]*1.2
