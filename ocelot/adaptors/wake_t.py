"""
This module contains methods for coupling the plasma simulation code Wake-T
with Ocelot

"""

import numpy as np

from ocelot.cpbd.beam import ParticleArray
from ocelot.common.globals import m_e_GeV


def waket_beam_to_parray(waket_beam, gamma_ref=None, z_ref=None):
    """
    Converts a Wake-T particle beam to an Ocelot ParticleArray.

    Parameters:
    -----------
    waket_beam : ParticleBunch (Wake-T class)

    gamma_ref : float

    z_ref : float


    Returns:
    --------

    """
    # Extract particle coordinates.
    x = waket_beam.x   # [m]
    y = waket_beam.y   # [m]
    z = waket_beam.xi  # [m]
    px = waket_beam.px # [m_e * c]
    py = waket_beam.py # [m_e * c]
    pz = waket_beam.pz # [m_e * c]
    q = waket_beam.q   # [C]

    # Calculate gamma.
    gamma = np.sqrt(1 + px**2 + py**2 + pz**2)

    # Determine reference energy and longitudinal position.
    if gamma_ref is None:
        gamma_ref = np.average(gamma, weights=q)
    if z_ref is None:
        z_ref = np.average(z)

    # Calculate momentum deviation (dp) and kinetic momentum (p_kin).
    b_ref = np.sqrt(1 - gamma_ref**(-2))
    dp = (gamma-gamma_ref)/(gamma_ref*b_ref)
    p_kin = np.sqrt(gamma**2 - 1)

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
    