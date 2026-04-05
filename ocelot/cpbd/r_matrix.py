__author__ = 'Sergey Tomin'

import logging
import numpy as np

from ocelot.common.globals import m_e_GeV, speed_of_light

logger = logging.getLogger(__name__)


def _broadcast_base(*values, xp=np):
    base = xp.asarray(0.0)
    for value in values:
        base = base + xp.asarray(value)
    return base


def _relativistic_terms(energy, xp=np):
    gamma = xp.asarray(energy) / m_e_GeV
    nonzero_gamma = gamma != 0
    gamma_safe = xp.where(nonzero_gamma, gamma, xp.ones_like(gamma))
    igamma2 = xp.where(nonzero_gamma, 1.0 / (gamma_safe * gamma_safe), xp.zeros_like(gamma))
    beta = xp.sqrt(1.0 - igamma2)
    return beta, igamma2


def _focusing_terms(z, k2, xp=np):
    z = xp.asarray(z)
    k2 = xp.asarray(k2)
    abs_k2 = xp.abs(k2)
    nonzero_k2 = k2 != 0
    positive_k2 = k2 >= 0
    k = xp.sqrt(abs_k2)
    k_safe = xp.where(nonzero_k2, k, xp.ones_like(k))
    k2_safe = xp.where(nonzero_k2, k2, xp.ones_like(k2))
    phase = k * z

    c = xp.where(positive_k2, xp.cos(phase), xp.cosh(phase))
    s = xp.where(positive_k2, xp.sin(phase) / k_safe, xp.sinh(phase) / k_safe)
    d = (1.0 - c) / k2_safe
    z_minus_s_over_k2 = (z - s) / k2_safe

    c = xp.where(nonzero_k2, c, xp.ones_like(c))
    s = xp.where(nonzero_k2, s, z)
    d = xp.where(nonzero_k2, d, 0.5 * z * z)
    z_minus_s_over_k2 = xp.where(nonzero_k2, z_minus_s_over_k2, z * z * z / 6.0)
    return c, s, d, z_minus_s_over_k2


def rot_mtx(angle, xp=np):
    if xp is np:
        cs = np.cos(angle)
        sn = np.sin(angle)
        return np.array([[cs, 0., sn, 0., 0., 0.],
                         [0., cs, 0., sn, 0., 0.],
                         [-sn, 0., cs, 0., 0., 0.],
                         [0., -sn, 0., cs, 0., 0.],
                         [0., 0., 0., 0., 1., 0.],
                         [0., 0., 0., 0., 0., 1.]])

    angle = xp.asarray(angle)
    cs = xp.cos(angle)
    sn = xp.sin(angle)
    zero = xp.zeros_like(cs)
    one = xp.ones_like(cs)
    return xp.stack([
        xp.stack([cs, zero, sn, zero, zero, zero]),
        xp.stack([zero, cs, zero, sn, zero, zero]),
        xp.stack([-sn, zero, cs, zero, zero, zero]),
        xp.stack([zero, -sn, zero, cs, zero, zero]),
        xp.stack([zero, zero, zero, zero, one, zero]),
        xp.stack([zero, zero, zero, zero, zero, one]),
    ])


def linear_magnet_matrix(z, k1, hx, sum_tilts=0.0, energy=0.0, xp=np):
    beta, igamma2 = _relativistic_terms(energy, xp=xp)
    beta2 = beta * beta

    kx2 = xp.asarray(k1) + xp.asarray(hx) * xp.asarray(hx)
    ky2 = -xp.asarray(k1)
    cx, sx, dx_scale, r56_scale = _focusing_terms(z, kx2, xp=xp)
    cy, sy, _, _ = _focusing_terms(z, ky2, xp=xp)
    dx = xp.asarray(hx) * dx_scale
    r56 = xp.asarray(hx) * xp.asarray(hx) * r56_scale / beta2 - xp.asarray(z) * igamma2 / beta2

    base = _broadcast_base(cx, cy, sx, sy, dx, r56, beta, xp=xp)
    zero = xp.zeros_like(base)
    one = xp.ones_like(base)
    u_matrix = xp.stack([
        xp.stack([cx, sx, zero, zero, zero, dx / beta]),
        xp.stack([-kx2 * sx, cx, zero, zero, zero, sx * xp.asarray(hx) / beta]),
        xp.stack([zero, zero, cy, sy, zero, zero]),
        xp.stack([zero, zero, -ky2 * sy, cy, zero, zero]),
        xp.stack([xp.asarray(hx) * sx / beta, dx / beta, zero, zero, one, r56]),
        xp.stack([zero, zero, zero, zero, zero, one]),
    ])

    if xp is np and sum_tilts == 0.0:
        return u_matrix

    r_tilt = rot_mtx(sum_tilts, xp=xp)
    return xp.matmul(rot_mtx(-sum_tilts, xp=xp), xp.matmul(u_matrix, r_tilt))


def bend_edge_matrix(h, edge, gap=0.0, fint=0.0, xp=np):
    h = xp.asarray(h)
    edge = xp.asarray(edge)
    gap = xp.asarray(gap)
    fint = xp.asarray(fint)
    sec_edge = 1.0 / xp.cos(edge)
    phi = fint * h * gap * sec_edge * (1.0 + xp.sin(edge) ** 2)

    base = _broadcast_base(h, edge, gap, fint, xp=xp)
    zero = xp.zeros_like(base)
    one = xp.ones_like(base)
    return xp.stack([
        xp.stack([one, zero, zero, zero, zero, zero]),
        xp.stack([h * xp.tan(edge), one, zero, zero, zero, zero]),
        xp.stack([zero, zero, one, zero, zero, zero]),
        xp.stack([zero, zero, -h * xp.tan(edge - phi), one, zero, zero]),
        xp.stack([zero, zero, zero, zero, one, zero]),
        xp.stack([zero, zero, zero, zero, zero, one]),
    ])


def cavity_coupler_edge_matrix(v, phi, vxx=0.0, vxy=0.0, energy=0.0, xp=np):
    phi_rad = xp.asarray(phi) * xp.pi / 180.0
    energy = xp.asarray(energy)
    valid_energy = energy != 0
    energy_safe = xp.where(valid_energy, energy, xp.ones_like(energy))
    rf_drive = xp.asarray(v) * xp.exp(1j * phi_rad)
    m21 = xp.real(xp.asarray(vxx) * rf_drive) / energy_safe
    m23 = xp.real(xp.asarray(vxy) * rf_drive) / energy_safe
    zero_kick = xp.zeros_like(m21 + m23)
    m21 = xp.where(valid_energy, m21, zero_kick)
    m23 = xp.where(valid_energy, m23, zero_kick)

    base = _broadcast_base(m21, m23, energy, xp=xp)
    zero = xp.zeros_like(xp.real(base))
    one = xp.ones_like(zero)
    return xp.stack([
        xp.stack([one, zero, zero, zero, zero, zero]),
        xp.stack([m21, one, m23, zero, zero, zero]),
        xp.stack([zero, zero, one, zero, zero, zero]),
        xp.stack([m23, zero, -m21, one, zero, zero]),
        xp.stack([zero, zero, zero, zero, one, zero]),
        xp.stack([zero, zero, zero, zero, zero, one]),
    ])


def standing_wave_cavity_matrix(z, voltage, energy, freq, phi=0.0, xp=np):
    z = xp.asarray(z)
    voltage = xp.asarray(voltage)
    energy = xp.asarray(energy)
    freq = xp.asarray(freq)
    phi_rad = xp.asarray(phi) * xp.pi / 180.0

    drift_matrix = uni_matrix(z, 0.0, hx=0.0, sum_tilts=0.0, energy=energy, xp=xp)

    de = voltage * xp.cos(phi_rad)
    gamma_i = energy / m_e_GeV
    gamma_f = (energy + de) / m_e_GeV
    valid_map = xp.logical_and(voltage != 0, xp.logical_and(energy != 0, gamma_f > 0))

    if xp is np and np.any(np.asarray((np.asarray(voltage) != 0) & (np.asarray(energy) == 0))):
        logger.error("CAVITY: Initial energy is 0, check ParticleArray.E or Twiss.E OR cavity.v must be 0")

    z_safe = xp.where(z != 0, z, xp.ones_like(z))
    gamma_i_safe = xp.where(valid_map, gamma_i, xp.ones_like(gamma_i))
    gamma_f_safe = xp.where(valid_map, gamma_f, xp.ones_like(gamma_f))
    delta_gamma = gamma_f_safe - gamma_i_safe
    delta_gamma_nonzero = xp.abs(delta_gamma) > 1.0e-10
    delta_gamma_safe = xp.where(delta_gamma_nonzero, delta_gamma, xp.ones_like(delta_gamma))

    cos_phi = xp.cos(phi_rad)
    sin_phi = xp.sin(phi_rad)
    cos_phi_nonzero = xp.abs(cos_phi) > 1.0e-14
    cos_phi_safe = xp.where(cos_phi_nonzero, cos_phi, xp.ones_like(cos_phi))
    log_ratio = xp.log(gamma_f_safe / gamma_i_safe)
    alpha = xp.where(cos_phi_nonzero, xp.sqrt(1.0 / 8.0) * log_ratio / cos_phi_safe, xp.zeros_like(log_ratio))
    sin_alpha = xp.sin(alpha)
    cos_alpha = xp.cos(alpha)

    ep = delta_gamma / z_safe
    ep_nonzero = xp.abs(ep) > 1.0e-10
    ep_safe = xp.where(ep_nonzero, ep, xp.ones_like(ep))

    r11 = cos_alpha - xp.sqrt(2.0) * cos_phi * sin_alpha
    r12 = xp.where(ep_nonzero, xp.sqrt(8.0) * gamma_i_safe / ep_safe * cos_phi * sin_alpha, z)
    r21 = -ep / gamma_f_safe * (cos_phi / xp.sqrt(2.0) + 1.0 / (xp.sqrt(8.0) * cos_phi_safe)) * sin_alpha
    r22 = gamma_i_safe / gamma_f_safe * (cos_alpha + xp.sqrt(2.0) * cos_phi * sin_alpha)

    beta_i_rel = xp.sqrt(xp.maximum(1.0 - 1.0 / (gamma_i_safe * gamma_i_safe), 0.0))
    beta_f_rel = xp.sqrt(xp.maximum(1.0 - 1.0 / (gamma_f_safe * gamma_f_safe), 0.0))
    beta_i_safe = xp.where(beta_i_rel != 0, beta_i_rel, xp.ones_like(beta_i_rel))
    beta_f_safe = xp.where(beta_f_rel != 0, beta_f_rel, xp.ones_like(beta_f_rel))

    k_rf = 2.0 * xp.pi * freq / speed_of_light
    r56_expr = -z / (gamma_f_safe * gamma_f_safe * gamma_i_safe * beta_f_safe) * (gamma_f_safe + gamma_i_safe) / (beta_f_safe + beta_i_safe)

    zero_energy_gain = xp.abs(delta_gamma) < 1.0e-8 * xp.abs(gamma_i_safe)
    is_zero_crossing = xp.logical_and(valid_map, xp.logical_and(zero_energy_gain, xp.abs(cos_phi) < 1.0e-3))
    general_den = beta_f_safe * gamma_f_safe * delta_gamma_safe * delta_gamma_safe
    r55_zero_expr = -k_rf * z * voltage / (2.0 * m_e_GeV * gamma_i_safe ** 3 * beta_i_safe ** 2)
    r55_general_expr = (
        k_rf * z * beta_i_safe * voltage / m_e_GeV * sin_phi
        * (gamma_i_safe * gamma_f_safe * (beta_i_safe * beta_f_safe - 1.0) + 1.0)
        / general_den
    )
    r55_cor = xp.where(valid_map, xp.where(is_zero_crossing, r55_zero_expr, r55_general_expr), xp.zeros_like(r55_general_expr))
    r65 = xp.where(valid_map, k_rf * sin_phi * voltage / (gamma_f_safe * beta_f_safe * m_e_GeV), xp.zeros_like(gamma_f_safe))
    r66 = xp.where(valid_map, gamma_i_safe / gamma_f_safe * beta_i_safe / beta_f_safe, xp.ones_like(gamma_f_safe))
    r56 = xp.where(valid_map, r56_expr, drift_matrix[4, 5])

    cavity_matrix = xp.stack([
        xp.stack([r11, r12, 0.0 * r11, 0.0 * r11, 0.0 * r11, 0.0 * r11]),
        xp.stack([r21, r22, 0.0 * r11, 0.0 * r11, 0.0 * r11, 0.0 * r11]),
        xp.stack([0.0 * r11, 0.0 * r11, r11, r12, 0.0 * r11, 0.0 * r11]),
        xp.stack([0.0 * r11, 0.0 * r11, r21, r22, 0.0 * r11, 0.0 * r11]),
        xp.stack([0.0 * r11, 0.0 * r11, 0.0 * r11, 0.0 * r11, 1.0 + r55_cor, r56]),
        xp.stack([0.0 * r11, 0.0 * r11, 0.0 * r11, 0.0 * r11, r65, r66]),
    ])
    return xp.where(valid_map, cavity_matrix, drift_matrix)


def uni_matrix(z, k1, hx, sum_tilts=0., energy=0., xp=np):
    """
    universal matrix. The function creates R-matrix from given parameters.
    r = element.l/element.angle
    +K - focusing lens, -K - defoc

    :param z: element length [m]
    :param k1: quadrupole strength [1/m**2]
    :param hx: the curvature (1/r) of the element [1/m]
    :param sum_tilts: rotation relative to longitudinal axis [rad]
    :param energy: the beam energy [GeV]
    :param xp: backend-compatible array namespace, e.g. numpy or jax.numpy
    :return: R-matrix [6, 6]
    """
    return linear_magnet_matrix(z, k1, hx, sum_tilts=sum_tilts, energy=energy, xp=xp)
