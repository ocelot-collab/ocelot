import beamphysics as pmd
import numpy as np

from ocelot.common.globals import m_e_eV
from ocelot.cpbd.beam import ParticleArray


def particle_group_to_parray(pgroup: pmd.ParticleGroup) -> ParticleArray:
    """Import an electron bunch with Cartesian momenta in eV/c.

    Choose the charge-weighted mean total energy as the reference. Derive
    electron energies from the input momenta with Ocelot's electron mass so
    that conversion back to Cartesian momenta preserves the mass shell.

    The longitudinal convention is tau = -z (metres). A ParticleGroup with
    arrival times at a screen must first be converted to a spatial bunch.
    This longitudinal convention assumes an ultrarelativistic bunch.
    """
    if pgroup.species != "electron":
        raise ValueError("The Ocelot ParticleArray adaptor requires electrons.")
    if np.any(pgroup.pz < 0):
        raise ValueError("Ocelot tracks forward-moving particles and requires pz >= 0.")
    parray = ParticleArray(pgroup.n_particle)
    energies = np.hypot(pgroup.p, m_e_eV)
    reference_energy = np.average(energies, weights=pgroup.weight)
    parray.E = reference_energy * 1e-9
    if parray.E <= m_e_eV * 1e-9:
        raise ValueError("The reference momentum must be positive.")
    reference_momentum = parray.p0c * 1e9

    rpart = parray.rparticles
    rpart[0] = pgroup.x
    rpart[1] = pgroup.px / reference_momentum
    rpart[2] = pgroup.y
    rpart[3] = pgroup.py / reference_momentum
    rpart[4] = -pgroup.z
    rpart[5] = (energies - reference_energy) / reference_momentum

    parray.q_array = pgroup.weight

    return parray


def load_pmd(filename: str) -> ParticleArray:
    return particle_group_to_parray(pmd.ParticleGroup(h5=filename))


def particle_array_to_particle_group(parray: ParticleArray) -> pmd.ParticleGroup:
    """Export Cartesian momenta in eV/c and a spatial bunch with z = -tau in metres."""
    px = parray.px() * parray.p0c * 1e9  # to eV
    py = parray.py() * parray.p0c * 1e9  # to eV
    pz = parray.pz * parray.p0c * 1e9  # to eV
    z = -parray.tau()

    data = {
        "x": parray.x(),
        "px": px,
        "y": parray.y(),
        "py": py,
        "z": z,
        "pz": pz,
        "t": np.zeros_like(px),
        "weight": parray.q_array,
        "status": np.ones_like(px),
        "species": "electron",
    }

    return pmd.ParticleGroup(data=data)
