from ocelot.cpbd.beam import ParticleArray
import numpy as np
import pytest

try:
    import beamphysics as pmd
except ImportError:
    IS_PMD_INSTALLED = False
else:
    IS_PMD_INSTALLED = True
    from ocelot.adaptors import pmd as pmd_adaptor

# Define decorator to skip tests if openpmd-beamphysics is not installed.
only_if_pmd_installed = pytest.mark.skipif(
    not IS_PMD_INSTALLED, reason='PMD required to run tests')


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
    # Compare physical momenta in eV/c. A comparison of normalized coordinates
    # alone would miss an inconsistent energy and momentum reference.
    np.testing.assert_allclose(pgroup.x, parray.x())
    np.testing.assert_allclose(pgroup.y, parray.y())
    np.testing.assert_allclose(pgroup.z, -parray.tau())
    np.testing.assert_allclose(pgroup.px, parray.px() * parray.p0c * 1e9)
    np.testing.assert_allclose(pgroup.py, parray.py() * parray.p0c * 1e9)
    np.testing.assert_allclose(pgroup.pz, parray.pz * parray.p0c * 1e9)

    np.testing.assert_allclose(pgroup.weight, parray.q_array)


@only_if_pmd_installed
def test_particle_group_to_parray():
    """Convertiong of ParticleGroup to ParticleArray"""
    # instantiate a ParticleGroup and make corresponding ParticleArray
    pgroup = pmd.ParticleGroup(data=PARTICLE_GROUP_DATA)
    parray = pmd_adaptor.particle_group_to_parray(pgroup)
    compare_particle_group_with_array(pgroup, parray)


@only_if_pmd_installed
def test_load_pmd(tmp_pmdh5):
    """Loading of PMD files"""
    parray = pmd_adaptor.load_pmd(str(tmp_pmdh5))
    pgroup = pmd.ParticleGroup(data=PARTICLE_GROUP_DATA)
    compare_particle_group_with_array(pgroup, parray)


@only_if_pmd_installed
def test_parray_to_particle_group(pmd_parray):
    """Conversion of ParticleArray to ParticleGroup"""
    pgroup = pmd_adaptor.particle_array_to_particle_group(pmd_parray)
    compare_particle_group_with_array(pgroup, pmd_parray)


@only_if_pmd_installed
def test_particle_array_round_trip_preserves_longitudinal_order():
    parray = ParticleArray(3)
    parray.E = 1.0
    parray.rparticles[:] = [
        [1e-5, -2e-5, 3e-5],
        [1e-6, 2e-6, -3e-6],
        [-4e-5, 5e-5, 6e-5],
        [-2e-6, 3e-6, 4e-6],
        [-3e-4, 2e-4, 7e-4],
        [-0.01, 0.0, 0.01],
    ]
    parray.q_array[:] = [1e-12, 2e-12, 1e-12]

    pgroup = pmd_adaptor.particle_array_to_particle_group(parray)
    restored = pmd_adaptor.particle_group_to_parray(pgroup)

    np.testing.assert_array_equal(restored.tau(), parray.tau())
    np.testing.assert_allclose(restored.rparticles, parray.rparticles, rtol=1e-9, atol=1e-14)
    np.testing.assert_array_equal(restored.q_array, parray.q_array)
    np.testing.assert_allclose(restored.E, parray.E, rtol=1e-12)


@only_if_pmd_installed
def test_low_energy_bunch_keeps_each_cartesian_momentum():
    data = dict(PARTICLE_GROUP_DATA, px=[1e4, 2e4, 0.], py=[-3e4, 4e4, 0.],
                pz=[1e5, 1e6, 2e6], weight=[1e-12, 2e-12, 1e-12])
    original = pmd.ParticleGroup(data=data)
    parray = pmd_adaptor.particle_group_to_parray(original)
    restored = pmd_adaptor.particle_array_to_particle_group(parray)
    for component in ("px", "py", "pz"):
        np.testing.assert_allclose(restored[component], original[component], rtol=3e-14)
    np.testing.assert_allclose(np.average(parray.p(), weights=parray.q_array), 0., atol=1e-15)


@only_if_pmd_installed
@pytest.mark.parametrize("changes, message", [
    ({"species": "proton"}, "electrons"),
    ({"pz": [1e9, -1e9, 1e9]}, "pz"),
    ({"pz": [0., 0., 0.], "px": [0., 0., 0.], "py": [0., 0., 0.]}, "reference momentum"),
])
def test_unrepresentable_group_rejected(changes, message):
    group = pmd.ParticleGroup(data=dict(PARTICLE_GROUP_DATA, **changes))
    with pytest.raises(ValueError, match=message):
        pmd_adaptor.particle_group_to_parray(group)
