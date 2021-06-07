from typing import Optional

import pmd_beamphysics as pmd
import numpy as np

from ocelot.cpbd.beam import ParticleArray


def particle_group_to_parray(pgroup: pmd.ParticleGroup) -> ParticleArray:
    """Construct an Ocelot ParticleArray from an openPMD-beamphysics ParticleGroup.
    The particle type is assumed to be electrons.

    :param pgroup: ParticleGroup from which to construct the ParticleArray
    :return: ParticleArray corresponding to the provided ParticleGroup
    :rtype: ParticleArray

    """
    parray = ParticleArray(pgroup.n_particle)

    # These are in eV, but so are px and py so this is OK.
    reference_momentum = pgroup.avg("p")
    reference_energy = pgroup.avg("energy")

    rpart = parray.rparticles
    rpart[0] = pgroup.x
    rpart[1] = pgroup.px / reference_momentum
    rpart[2] = pgroup.y
    rpart[3] = pgroup.py / reference_momentum
    rpart[4] = pgroup.z
    rpart[5] = (pgroup.energy - reference_energy) / reference_momentum

    parray.E = reference_energy * 1e-9  # Convert eV to GeV
    parray.q_array = pgroup.weight

    return parray


def load_pmd(filename: str) -> ParticleArray:
    return particle_group_to_parray(pmd.ParticleGroup(h5=filename))


def particle_array_to_particle_group(parray: ParticleArray) -> pmd.ParticleGroup:
    px = parray.px() * parray.p0 * 1e9  # to eV
    py = parray.py() * parray.p0 * 1e9  # to eV
    pz = parray.pz * parray.p0 * 1e9  # to eV

    data = {
        "x": parray.x(),
        "px": px,
        "y": parray.y(),
        "py": py,
        "z": parray.tau(),
        "pz": pz,
        "t": np.zeros_like(px),
        "weight": parray.q_array,
        "status": np.ones_like(px),
        "species": "electron",
    }

    return pmd.ParticleGroup(data=data)
