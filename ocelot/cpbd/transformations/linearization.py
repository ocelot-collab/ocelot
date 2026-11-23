"""Local derivatives of active single-particle tracking maps."""

import numpy as np

from ocelot.common.globals import m_e_GeV, speed_of_light
from ocelot.cpbd.rf_maps import rf_long_t_terms
from ocelot.cpbd.transformations import CavityTM, SecondTM, TMTypes, TransferMap, TWCavityTM
from ocelot.cpbd.transformations.transformation import Transformation


def map_derivatives(
    tm: Transformation,
    energy: float,
    orbit: np.ndarray,
    dispersion: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Return the Jacobian and, when requested, its derivative along dispersion."""
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
        probes = np.column_stack((point[:, None] + offsets, point[:, None] - offsets))
        mapped = tm.map_function(probes, energy)
        return (mapped[:, :6] - mapped[:, 6:]) / (2. * step)

    J = jacobian(orbit)
    if dispersion is None:
        return J, None
    norm = np.max(np.abs(dispersion))
    if norm == 0.:
        return J, np.zeros_like(J)
    delta = step / norm
    shift = delta * dispersion
    return J, (jacobian(orbit + shift) - jacobian(orbit - shift)) / (2. * delta)
