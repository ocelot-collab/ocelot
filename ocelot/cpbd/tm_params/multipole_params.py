import numpy as np

from ocelot.common.globals import m_e_GeV
from ocelot.cpbd.tm_params.tm_params import TMParams


def thin_dipole_dispersion(angle: complex, energy: float) -> complex:
    """Thin-bend angular dispersion for delta = (E-E0)/(p0*c)."""
    beta = np.sqrt(1. - (m_e_GeV / energy)**2) if energy else 1.
    return angle / beta


def thin_multipole_matrix(dipole: complex, gradient: complex) -> np.ndarray:
    """Thin-kick Jacobian in magnet axes for the local complex field gradient."""
    R = np.eye(6)
    R[1, 0] = -gradient.real
    R[3, 2] = gradient.real
    R[1, 2] = R[3, 0] = gradient.imag
    # The conjugate momentum of arrival time tau is -delta.
    R[1, 5] = R[4, 0] = dipole.real
    R[3, 5] = R[4, 2] = -dipole.imag
    return R


class MultipoleParams(TMParams):
    """
    Coefficient list for the dedicated ``MultipoleTM`` kick polynomial.

    ``kn`` stores the normal multipole strengths in the order expected by the
    multipole transformation.
    """

    def __init__(self, kn, dx=0., dy=0., tilt=0.) -> None:
        super().__init__()
        self.kn = kn
        self.dx = dx
        self.dy = dy
        self.tilt = tilt
