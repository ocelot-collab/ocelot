__author__ = 'Sergey Tomin'
"""
module contains lat2input function which creates python input string
(ocelot lattice) for a lattice object
author sergey.tomin
"""
import os
import numpy as np

from ocelot.adaptors.astra2ocelot import astraBeam2particleArray, particleArray2astraBeam
from ocelot.adaptors.csrtrack2ocelot import csrtrackBeam2particleArray, particleArray2csrtrackBeam

try:
    from ocelot.adaptors.pmd import load_pmd, particle_array_to_particle_group
except ImportError:
    pass

from ocelot.cpbd.beam import ParticleArray, Twiss, Beam


def save_particle_array2npz(filename, p_array):
    np.savez_compressed(filename, rparticles=p_array.rparticles,
                        q_array=p_array.q_array,
                        E=p_array.E, s=p_array.s)


def load_particle_array_from_npz(filename, print_params=False):
    """
    Load beam file in npz format and return ParticleArray

    :param filename:
    :return:
    """
    p_array = ParticleArray()
    with np.load(filename) as data:
        for key in data.keys():
            p_array.__dict__[key] = data[key]
    return p_array


def load_particle_array(filename, print_params=False):
    """
    Universal function to load beam file, *.ast (ASTRA), *.fmt1 (CSRTrack) or *.npz format

    Note that downloading ParticleArray from the astra file (.ast) and saving it back does not give the same distribution.
    The difference arises because the array of particles does not have a reference particle, and in this case
    the first particle is used as a reference.

    :param filename: path to file, filename.ast or filename.npz
    :return: ParticleArray
    """
    name, file_extension = os.path.splitext(filename)
    if file_extension == ".npz":
        parray = load_particle_array_from_npz(filename, print_params=False)
    elif file_extension in [".ast", ".001"]:
        parray = astraBeam2particleArray(filename, print_params=False)
    elif file_extension in [".fmt1"]:
        parray = csrtrackBeam2particleArray(filename)
    elif file_extension == ".h5":
        parray = load_pmd(filename)
    else:
        raise Exception("Unknown format of the beam file: " + file_extension + " but must be *.ast, *fmt1 or *.npz ")

    if print_params:
        print(parray)

    return parray


def save_particle_array(filename, p_array):
    """
    Universal function to save beam file, *.ast (ASTRA), *.fmt1 (CSRTrack) or *.npz format

    Note that downloading ParticleArray from the astra file (.ast) and saving it back does not give the same distribution.
    The difference arises because the array of particles does not have a reference particle, and in this case
    the first particle is used as a reference.

    :param filename: path to file, filename.ast or filename.npz
    :return: ParticleArray
    """
    name, file_extension = os.path.splitext(filename)
    if file_extension == ".npz":
        save_particle_array2npz(filename, p_array)
    elif file_extension == ".ast":
        particleArray2astraBeam(p_array, filename)
    elif file_extension == ".fmt1":
        particleArray2csrtrackBeam(p_array, filename)
    elif file_extension == ".h5":
        particle_array_to_particle_group(p_array).write(filename)
    else:
        particle_array_to_particle_group(p_array).write(filename)
        raise Exception("Unknown format of the beam file: " + file_extension + " but must be *.ast, *.fmt1 or *.npz")

