"""
This module contains methods for coupling the plasma simulation code Wake-T
with Ocelot.

"""

import numpy as np
try:
    from wake_t.driver_witness import ParticleBunch
    wake_t_installed = True
except ImportError:
    wake_t_installed = False

from ocelot.cpbd.beam import ParticleArray
from ocelot.common.globals import m_e_GeV


def wake_t_beam_to_parray(wake_t_beam, gamma_ref=None, z_ref=None):
    """
    Converts a Wake-T particle beam to an Ocelot ParticleArray.

    Parameters:
    -----------
    wake_t_beam : ParticleBunch (Wake-T class)
        The original particle distribution from Wake-T.

    gamma_ref : float (Optional)
        Reference energy of the particle beam used for tracking in Ocelot. If
        not specified, the reference energy will be taken as the average
        energy of the input distribution.

    z_ref : float (Optional)
        Reference longitudinal position of the particle beam used for tracking
        in Ocelot. If not specified, the reference value will be taken as the
        average longitudinal position of the input distribution.

    Returns:
    --------
    An Ocelot ParticleArray.

    """
    # Extract particle coordinates.
    x = wake_t_beam.x  # [m]
    y = wake_t_beam.y  # [m]
    z = wake_t_beam.xi + wake_t_beam.prop_distance  # [m]
    px = wake_t_beam.px  # [m_e * c]
    py = wake_t_beam.py  # [m_e * c]
    pz = wake_t_beam.pz  # [m_e * c]
    q = wake_t_beam.q  # [C]

    # Calculate gamma.
    gamma = np.sqrt(1 + px**2 + py**2 + pz**2)

    # Determine reference energy and longitudinal position.
    if gamma_ref is None:
        gamma_ref = np.average(gamma, weights=q)
    if z_ref is None:
        z_ref = wake_t_beam.prop_distance

    # Calculate momentum deviation (dp) and kinetic momentum (p_kin).
    b_ref = np.sqrt(1 - gamma_ref**(-2))
    dp = (gamma-gamma_ref)/(gamma_ref*b_ref)
    p_kin = np.sqrt(gamma_ref**2 - 1)

    # Create particle array
    p_array = ParticleArray(len(q))
    p_array.rparticles[0] = x
    p_array.rparticles[2] = y
    p_array.rparticles[4] = -(z - z_ref)
    p_array.rparticles[1] = px/p_kin
    p_array.rparticles[3] = py/p_kin
    p_array.rparticles[5] = dp
    p_array.q_array[:] = q
    p_array.s = z_ref
    p_array.E = gamma_ref * m_e_GeV

    return p_array


def parray_to_wake_t_beam(p_array):
    """
    Converts an Ocelot ParticleArray to a Wake-T ParticleBunch.

    Parameters:
    -----------
    p_array : ParticleArray
        The Ocelot distribution to be converted.

    Returns:
    --------
    A Wake-T ParticleBunch.

    """
    # Check if Wake-T is installed.
    if not wake_t_installed:
        raise ImportError('Wake-T is not installed. '
                          'Cannot perform conversion to Wake-T ParticleBunch.')

    # Get beam matrix and reference values.
    beam_matrix = p_array.rparticles
    gamma_ref = p_array.E / m_e_GeV
    z_ref = p_array.s

    # Calculate gamma and kinetic momentum (p_kin).
    dp = beam_matrix[5]
    b_ref = np.sqrt(1 - gamma_ref**(-2))
    gamma = dp*gamma_ref*b_ref + gamma_ref
    p_kin = np.sqrt(gamma_ref**2 - 1)

    # Create coordinate arrays in Wake-T units.
    x = beam_matrix[0]
    px = beam_matrix[1] * p_kin
    y = beam_matrix[2]
    py = beam_matrix[3] * p_kin
    z = - beam_matrix[4]
    pz = np.sqrt(gamma**2 - px**2 - py**2 - 1)
    q = p_array.q_array

    # Create and return Wake-T distribution.
    return ParticleBunch(q, x, y, z, px, py, pz, prop_distance=z_ref)
