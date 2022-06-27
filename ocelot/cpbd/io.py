__author__ = 'Sergey Tomin'
"""
module contains lat2input function which creates python input string
(ocelot lattice) for a lattice object
author sergey.tomin
"""
import os
import logging
from typing import Union, List, Tuple, Optional, Any, Iterable, Type
from types import TracebackType
from pathlib import Path

import numpy as np


from ocelot.adaptors.astra2ocelot import (astraBeam2particleArray,
                                          particleArray2astraBeam)
from ocelot.adaptors.csrtrack2ocelot import (csrtrackBeam2particleArray,
                                             particleArray2csrtrackBeam)

try:
    from ocelot.adaptors.pmd import load_pmd, particle_array_to_particle_group
except ImportError:
    pass

from ocelot.cpbd.beam import ParticleArray, Twiss, Beam

try:
    from mpi4py import MPI
except ImportError:
    pass

try:
    import h5py
except ImportError:
    pass

_logger = logging.getLogger(__name__)

def is_an_mpi_process():
    return "OMPI_COMM_WORLD_SIZE" in os.environ

class HDF5FileWithMaybeMPI(h5py.File):
    def __init__ (self, *args, **kwargs):
        if is_an_mpi_process():
            _logger.debug("Opening h5py with mpi enabled")
            kwargs.update({"driver": "mpio", "comm": MPI.COMM_WORLD})
        # self.f = h5py.File(*args, **kwargs)
        super().__init__(*args, **kwargs)


class ParameterScanFile:
    def __init__(self, *args, **kwargs):
        self._hdf = HDF5FileWithMaybeMPI(*args, **kwargs)

    def __enter__(self):
        return self

    def __exit__(self, exctype: Optional[Type[BaseException]],
                 excinst: Optional[BaseException],
                 exctb: Optional[TracebackType]) -> bool:
        self.close()

    def close(self):
        self._hdf.close()

    @property
    def filename(self):
        return self._hdf.filename

    def new_run_output(self, pval, parray0, marker_names):
        next_run_name = self.next_run_name()
        self._initialise_one_root_level_group(next_run_name, pval, parray0, marker_names)
        return next_run_name

    def _write_parray(self, scan_index, parray_group_name, parray):
        run_name = self.run_names[scan_index]
        group = self._hdf[run_name][parray_group_name]
        self._write_beam_to_group(group, parray)

    def write_parray0(self, scan_index, parray0):
        self._write_parray(scan_index, "parray0", parray0)

    def write_parray1(self, scan_index, parray1):
        self._write_parray(scan_index, "parray1", parray1)

    def write_parray_marker(self, scan_index, marker_name, parray):
        self._write_parray(scan_index, f"markers/{marker_name}", parray)

    def _yield_parrays_from_run_group(self, from_run_address):
        # from_run address is the address relative to the groups named run-0,
        # run-1, etc.  so "parray0" for example not "run-1/parray0"!
        for run_name in self.run_names:
            yield h5py_group_to_parray(self._hdf[f"{run_name}/{from_run_address}"])

    def parray0s(self):
        yield from self._yield_parrays_from_run_group("parray0")

    def parray1s(self):
        yield from self._yield_parrays_from_run_group("parray1")

    def parray_markers(self, marker_name):
        if marker_name not in self.marker_names:
            raise ValueError(f"Marker name not found in file: {marker_name}")
        yield from self._yield_parrays_from_run_group(f"markers/{marker_name}")

    @property
    def run_names(self):
        keys = self._hdf.keys()
        return [key for key in keys if key.startswith("run-")]

    @property
    def marker_names(self):
        result = set()
        for run_name in self.run_names:
            result.update(list(self._hdf[run_name]["markers"]))
        return result

    @property
    def parameter_name(self):
        return self._hdf["parameter-name"][()].decode("utf-8")

    @parameter_name.setter
    def parameter_name(self, value):
        self._hdf["parameter-name"] = np.string_(value)

    @property
    def parameter_values(self):
        result = []
        for run_name in self.run_names:
            result.append(self._hdf[f"{run_name}/parameter-value"][()])
        return result

    def next_run_name(self):
        run_numbers = [int(n.split("-")[1]) for n in self.run_names]
        if not run_numbers:
            return "run-0"
        return self._run_index_to_name(max(run_numbers) + 1)

    def init_from_parameter_scanner(self, parameter_scanner) -> None:
        parameter_values = parameter_scanner.parameter_values
        parray0s = parameter_scanner.parray0
        if isinstance(parray0s, ParticleArray):
            parray0s = len(parameter_values) * [parray0s]
        marker_names = [marker.id for marker in parameter_scanner.markers]
        # parameter_values = pscanner.parameter_values
        # You need to make EVERY group and EVERY dataset in EVERY core when using
        # MPI-enabled h5py.  This is documented but also if you don't do it, it
        # won't fail loudly, you'll just get weird behaviour (e.g. all datasets the
        # same, blank, etc..)

        self.parameter_name = parameter_scanner.parameter_name

        assert not isinstance(parray0s, ParticleArray)

        for parameter, parray0 in zip(parameter_values, parray0s):
            self.new_run_output(parameter, parray0, marker_names)

    def _initialise_one_root_level_group(self, run_name: str,
                                         parameter: Any,
                                         parray0: ParticleArray,
                                         marker_names) -> None:
        """Initialises a group of the kind found at the root level of the file and
        everything within it."""
        group = self._hdf.create_group(run_name)
        g0 = group.create_group("parray0")
         # init with specific parray0 in case they each have different lengths.
        self._init_h5_parray_group(g0, parray0)
        g1 = group.create_group("parray1")
        self._init_h5_parray_group(g1, parray0)
        group["parameter-value"] = parameter

        # Initialise marker parray output.
        gm = group.create_group("markers")
        for marker_name in marker_names:
            self._init_h5_parray_group(gm.create_group(marker_name), parray0)

    @staticmethod
    def _run_index_to_name(run_index: int) -> str:
        return f"run-{run_index}"

    def _write_output(self, hdf, run_name: str, parray0: ParticleArray, parray1:
                      ParticleArray) -> None:
        # run_index provides where the particle array will be written.
        self._write_beam_to_group(hdf[f"{run_name}/parray0"], parray0)
        self._write_beam_to_group(hdf[f"{run_name}/parray1"], parray1)

    @staticmethod
    def _init_h5_parray_group(group, parray0: ParticleArray) -> None:
        # Have to explicitly set the dtypes.
        group.create_dataset("rparticles",
                             shape=parray0.rparticles.shape,
                             dtype=parray0.rparticles.dtype)
        group.create_dataset("q_array", shape=parray0.q_array.shape,
                             dtype=parray0.q_array.dtype)
        group.create_dataset("E", shape=(), dtype=type(parray0.E))
        group.create_dataset("s", shape=(), dtype=type(parray0.s))

    @staticmethod
    def _write_beam_to_group(group, parray: ParticleArray) -> None:
        group["rparticles"][:] = parray.rparticles
        group["q_array"][:] = parray.q_array
        group["E"][()] = parray.E
        group["s"][()] = parray.s


def h5py_group_to_parray(group):
    parray = ParticleArray()
    for key in 'E', 'q_array', 'rparticles', 's':
        setattr(parray, key, group[key][()])
    return parray

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
    p_array.lost_particle_recorder = p_array.LostParticleRecorder(p_array.n)
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
