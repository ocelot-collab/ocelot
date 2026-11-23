"""Local derivatives of active single-particle tracking maps."""

from math import factorial

import numpy as np

from ocelot.common.globals import m_e_GeV, speed_of_light
from ocelot.cpbd.r_matrix import rot_mtx
from ocelot.cpbd.rf_maps import rf_long_t_terms
from ocelot.cpbd.tm_params.multipole_params import thin_dipole_dispersion, thin_multipole_matrix
from ocelot.cpbd.transformations import CavityTM, MultipoleTM, SecondTM, TMTypes, TransferMap, TWCavityTM
from ocelot.cpbd.transformations.transformation import Transformation


def map_derivatives(
    tm: Transformation,
    energy: float,
    orbit: np.ndarray,
    dispersion: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Return the Jacobian and, when requested, its derivative along dispersion."""
    if isinstance(tm, MultipoleTM):
        params = tm.get_params(energy)
        rotation = rot_mtx(params.tilt)
        offset = np.array([params.dx, 0., params.dy, 0., 0., 0.])
        local = rotation @ (orbit - offset)
        z = local[0] + 1j * local[2]
        # Differentiate the thin kick polynomial in z = x + iy. A finite
        # difference would give even a centered octupole artificial focusing.
        gradient = sum(k * z**n / factorial(n) for n, k in enumerate(params.kn[1:]))
        dipole = thin_dipole_dispersion(params.kn[0], energy)
        J = rotation.T @ thin_multipole_matrix(dipole, gradient) @ rotation
        if dispersion is None:
            return J, None
        local_dispersion = rotation @ dispersion
        curvature = sum(k * z**n / factorial(n) for n, k in enumerate(params.kn[2:]))
        change = curvature * (local_dispersion[0] + 1j * local_dispersion[2])
        return J, rotation.T @ (thin_multipole_matrix(0., change) - np.eye(6)) @ rotation
    if isinstance(tm, SecondTM):
        params = tm.get_params(energy)
        tensor = params.get_rotated_T() if params.tilt else params.T
        hessian = tensor + tensor.swapaxes(1, 2)
        dJ = hessian @ dispersion if dispersion is not None else None
        return params.get_rotated_R() + hessian @ orbit, dJ
    if isinstance(tm, (CavityTM, TWCavityTM)):
        params = tm.get_params(energy)
        J = params.get_rotated_R()
        dJ = np.zeros_like(J) if dispersion is not None else None
        if tm.tm_type == TMTypes.MAIN and (energy or params.v):
            z = tm.length if tm.delta_length is None else tm.delta_length
            voltage = params.v * z / tm.length if tm.length else params.v
            phi = np.deg2rad(params.phi)
            t566, t556, t555 = rf_long_t_terms(z, voltage, energy, params.freq, phi)
            h_tau = np.array([[2. * t555, t556], [t556, 2. * t566]])
            J[4, 4:] += h_tau @ orbit[4:]
            # rf_longitudinal_kick uses the full sinusoid for outgoing delta.
            k = 2. * np.pi * params.freq / speed_of_light
            phase = phi - k * orbit[4]
            e_out = energy + voltage * np.cos(phi)
            pc_out = np.sqrt((e_out - m_e_GeV) * (e_out + m_e_GeV))
            J[5, 4] = voltage * k * np.sin(phase) / pc_out
            if dispersion is not None:
                dJ[4, 4:] = h_tau @ dispersion[4:]
                dJ[5, 4] = -voltage * k**2 * np.cos(phase) / pc_out * dispersion[4]
        return J, dJ
    if type(tm) is TransferMap:
        dJ = np.zeros((6, 6)) if dispersion is not None else None
        return tm.get_params(energy).get_rotated_R(), dJ

    # Non-polynomial tracking maps need local derivatives of their actual map
    # function. Evaluate coordinate probes together, without changing magnets
    # or the design reference energy. Scale the outer step in phase space.
    step = np.cbrt(np.finfo(float).eps)
    offsets = step * np.eye(6)

    def jacobian(point: np.ndarray) -> np.ndarray:
        # A two-point derivative mistakes a cubic kick for a quadrupole on
        # the octupole axis. Cancel that term with the five-point stencil.
        probes = point[:, None] + np.column_stack((offsets, -offsets, 2 * offsets, -2 * offsets))
        mapped = tm.map_function(probes, energy)
        plus, minus, far_plus, far_minus = np.split(mapped, 4, axis=1)
        return (8 * (plus - minus) - (far_plus - far_minus)) / (12 * step)

    J = jacobian(orbit)
    if dispersion is None:
        return J, None
    norm = np.max(np.abs(dispersion))
    if norm == 0.:
        return J, np.zeros_like(J)
    delta = step / norm
    shift = delta * dispersion
    return J, (jacobian(orbit + shift) - jacobian(orbit - shift)) / (2. * delta)
