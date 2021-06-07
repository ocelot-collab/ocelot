import numpy as np
import pmd_beamphysics as pmd
import pytest

from ocelot.adaptors import pmd as pmd_adaptor
from ocelot.common.globals import m_e_eV
from ocelot.cpbd.beam import ParticleArray


# Toy data for instantiating a ParticleGroup
PARTICLE_GROUP_DATA = {
    "x": np.array([1e-6, 2e-6, 3e-6]),
    "px": np.array([10, 20, 30]),
    "y": np.array([1e-7, 2e-7, 3e-7]),
    "py": np.array([15, 25, 35]),
    "z": np.array([0.1, 0.2, 0.3]),
    "pz": np.array([1e9, 1.01e9, 1.02e9]),  # i.e. ~1 GeV/c
    "t": np.array([1e-9, 2e-9, 3e-9]),
    "weight": np.array([1e-15, 1e-15, 1e-15]),
    "status": np.array([1, 1, 1]),
    "species": "electron",
}


@pytest.fixture
def tmp_pmdh5(tmp_path):
    """tmp path to written pmd h5 file defined from above PARTICLE_GROUP_DATA"""
    pg = pmd.ParticleGroup(data=PARTICLE_GROUP_DATA)
    pmd_path = tmp_path / "pmd.h5"
    pg.write(str(pmd_path))
    yield pmd_path


@pytest.fixture
def pmd_parray(tmp_pmdh5):
    """Ocelot ParticleArray fixture from the above PARTICLE_GROUP_DATA."""
    pgroup = pmd.ParticleGroup(h5=str(tmp_pmdh5))
    yield pmd_adaptor.particle_group_to_parray(pgroup)


def compare_particle_group_with_array(pgroup, parray):
    # Use the values from the pgroup for consistency as these are what
    # particle_group_to_parray uses internally.
    refmom = pgroup.avg("p")
    refenergy = pgroup.avg("energy")

    np.testing.assert_allclose(pgroup.x, parray.x())
    np.testing.assert_allclose(pgroup.px / refmom, parray.px()),
    np.testing.assert_allclose(pgroup.y, parray.y())
    np.testing.assert_allclose(pgroup.py / refmom, parray.py())

    np.testing.assert_allclose(pgroup.z, parray.tau())
    dp = (pgroup.energy - refenergy) / refmom
    np.testing.assert_allclose(dp, parray.p())

    np.testing.assert_allclose(pgroup.weight, parray.q_array)


def test_particle_group_to_parray():
    """Convertiong of ParticleGroup to ParticleArray"""
    # instantiate a ParticleGroup and make corresponding ParticleArray
    pgroup = pmd.ParticleGroup(data=PARTICLE_GROUP_DATA)
    parray = pmd_adaptor.particle_group_to_parray(pgroup)
    compare_particle_group_with_array(pgroup, parray)


def test_load_pmd(tmp_pmdh5):
    """Loading of PMD files"""
    parray = pmd_adaptor.load_pmd(str(tmp_pmdh5))
    pgroup = pmd.ParticleGroup(data=PARTICLE_GROUP_DATA)
    compare_particle_group_with_array(pgroup, parray)


def test_parray_to_particle_group(pmd_parray):
    """Conversion of ParticleArray to ParticleGroup"""
    pgroup = pmd_adaptor.particle_array_to_particle_group(pmd_parray)
    compare_particle_group_with_array(pgroup, pmd_parray)
