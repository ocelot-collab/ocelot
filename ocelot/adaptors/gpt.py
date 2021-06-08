from pathlib import Path
from typing import Mapping, Sequence, Tuple, Union

from ocelot.cpbd.beam import ParticleArray
from ocelot.common.globals import m_e_GeV, speed_of_light

import numpy as np
import easygdf

ScreenToutType = Mapping[str, Union[float, Sequence[float]]]


def touts_and_screens_to_particle_arrays(
    filename: str,
) -> Tuple[Tuple[float, ParticleArray], Tuple[float, ParticleArray]]:
    """Convert GDF file's touts and screens to postions/times with ParticleArray
    instances.

    :param filename: The path of the GDF file.

    """

    g = easygdf.load_screens_touts(filename)
    touts = g["touts"]
    touts = [(tout["time"], tout_to_particle_array(tout)) for tout in touts]

    screens = g["screens"]
    screens = [(scrn["position"], screen_to_particle_array(scrn)) for scrn in screens]

    return touts, screens


def tout_to_particle_array(tout: ScreenToutType) -> ParticleArray:
    """Convert a map of tout keys/value pairs to a ParticleArray.

    :param tout: The Mapping of tout variable names to values.

    """
    # This is an assumption, that z in the tout/screen coordinate system is the same as
    # z in the familiar accelerator physics coordinate system.
    dz = tout["z"] - np.mean(tout["z"])
    # particles in front of the reference (positive z) corresponds to a negative time
    # and therefore a negative tau = ct
    tau = -dz
    return _common_implementation(tout, tau)


def screen_to_particle_array(screen: ScreenToutType) -> ParticleArray:
    """Convert a map of screen keys/value pairs to a ParticleArray.

    :param screen: The Mapping of screen variable names to values.

    """
    dt = np.mean(screen["t"]) - screen["t"]  # Negative if in front of the reference
    tau = speed_of_light * dt
    return _common_implementation(screen, tau)


def _common_implementation(screen_or_tout: ScreenToutType, tau: float) -> ParticleArray:
    """Convert screen or tout to a ParticleArray with the given tau coordinate"""
    sot = screen_or_tout

    energy = sot["G"] * m_e_GeV
    momentum = (energy ** 2 - m_e_GeV ** 2) ** 0.5
    reference_energy = np.mean(energy)
    reference_momentum = np.mean(momentum)

    parray = ParticleArray(len(energy))
    rpart = parray.rparticles
    rpart[0] = np.mean(sot["x"]) - sot["x"]
    rpart[1] = sot["Bx"] * momentum / reference_momentum
    rpart[2] = np.mean(sot["y"]) - sot["y"]
    rpart[3] = sot["By"] * momentum / reference_momentum
    rpart[4] = tau
    rpart[5] = (energy - reference_energy) / reference_momentum
    parray.q_array = sot["nmacro"] * abs(sot["q"])
    parray.E = reference_energy

    return parray
