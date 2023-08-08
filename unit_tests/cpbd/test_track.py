import copy

import pytest
import numpy as np

from ocelot.cpbd.io import ParameterScanFile
from ocelot.cpbd.track import ParameterScanner
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.elements import SBend, Marker
from ocelot.cpbd.transformations import SecondTM
from ocelot.cpbd.beam import ParticleArray

def assert_parray_equality(first, second):
    np.testing.assert_equal(first.E, second.E)
    np.testing.assert_equal(first.s, second.s)
    np.testing.assert_equal(first.q_array, second.q_array)
    np.testing.assert_equal(first.rparticles, second.rparticles)

PARAMETER_NAME = "parameter-name"
MARKER_NAMES = ["pytest-marker-name-1", "pytest-marker-name-2"]
PARAMETER_VALUES = [0.1, 0.3]
FILENAME = "test.hdf5"

try:
    import h5py
except ImportError:
    pytest.skip(reason="No h5py installed", allow_module_level=True)

@pytest.fixture
def navigator():
    b1 = SBend(l=0.5, angle=-0.5)
    cell = [Marker(MARKER_NAMES[0]), b1, Marker(MARKER_NAMES[1])]
    magnetic_lattice = MagneticLattice(cell, method={"global": SecondTM})
    navi = Navigator(magnetic_lattice)

    return navi

@pytest.fixture()
def pscanner(navigator):
    parray0s = [ParticleArray.random(15), ParticleArray.random(10)]

    markers = [ele for ele in navigator.lat.sequence if isinstance(ele, Marker)]

    pscanner = ParameterScanner(navigator,
                             PARAMETER_VALUES,
                             parray0s,
                             parameter_name=PARAMETER_NAME,
                             markers=markers)

    return pscanner

def test_ParameterScanner_init(navigator):
    parray0s = [ParticleArray.random(15), ParticleArray.random(10)]

    markers = [ele for ele in navigator.lat.sequence if isinstance(ele, Marker)]

    pscanner = ParameterScanner(navigator,
                             PARAMETER_VALUES,
                             parray0s,
                             parameter_name=PARAMETER_NAME,
                             markers=markers)

    assert pscanner.navigator == navigator
    assert pscanner.parameter_values == PARAMETER_VALUES
    assert pscanner.parray0 == parray0s
    assert pscanner.parameter_name == PARAMETER_NAME
    assert pscanner.markers == markers

    # initialised with only a single parray + default init.:
    pscanner = ParameterScanner(navigator,
                             PARAMETER_VALUES,
                             parray0s[0])

    assert pscanner.markers == []
    assert pscanner.parameter_name == "parameter"


def test_throws_when_bad_length_mismatch_init(navigator):
    parray0s = [ParticleArray.random(15), ParticleArray.random(10)]

    paramvals = (len(parray0s) + 1) * [0.3]  # Always longer than the number of parray0s

    with pytest.raises(ValueError):
        pscanner = ParameterScanner(navigator, paramvals, parray0s)

def test_njobs(pscanner):
    assert pscanner.njobs == 2

def test_ParameterScanner_prepare_navigator(pscanner, monkeypatch):
    # Assert its not equal, set deepcopy to simply not do a copy temporarily,
    # and then assert equality. this more or less guarantees prepare_navigator
    # works relatively easily.
    assert pscanner.prepare_navigator() != pscanner.navigator
    def nocopy(arg):
        return arg
    monkeypatch.setattr(copy, "deepcopy", nocopy)
    assert pscanner.prepare_navigator() == pscanner.navigator


def test_ParameterScanner_generate_unique_navigators_with_parray0s(pscanner):
    navis, parray0s = pscanner.generate_unique_navigators_with_parray0s()

    # Check that CopyBeam has been attached
    proclist0 = navis[0].process_table.proc_list
    proclist1 = navis[0].process_table.proc_list

    # One for each marker (there are two markers)
    assert len(proclist0) == 2
    assert len(proclist1) == 2

    assert proclist0[0].name == MARKER_NAMES[0]
    assert proclist0[1].name == MARKER_NAMES[1]
    assert proclist1[0].name == MARKER_NAMES[0]
    assert proclist1[1].name == MARKER_NAMES[1]

    assert_parray_equality(parray0s[0], pscanner.parray0[0])
    assert_parray_equality(parray0s[1], pscanner.parray0[1])

def test_ParameterScanner_prepare_track_payloads(pscanner):
    navis, parray0s = pscanner.generate_unique_navigators_with_parray0s()
    payloads = pscanner.prepare_track_payloads(navis, parray0s)

    assert navis[0].lat is payloads[0].lattice
    assert navis[1].lat is payloads[1].lattice

    assert_parray_equality(payloads[0].parray0, parray0s[0])
    assert_parray_equality(payloads[1].parray0, parray0s[1])

    assert navis[0] is payloads[0].navigator
    assert navis[1] is payloads[1].navigator

    assert payloads[0].job_index == 0
    assert payloads[1].job_index == 1

    assert payloads[0].calc_tws == False
    assert payloads[1].calc_tws == False

    assert payloads[0].bounds == None
    assert payloads[1].bounds == None

    assert payloads[0].return_df == True
    assert payloads[1].return_df == True

    assert payloads[0].overwrite_progress == False
    assert payloads[1].overwrite_progress == False

    assert payloads[0].kwargs == None
    assert payloads[1].kwargs == None

def test_ParameterScanner_scan(pscanner, tmp_path):
    pscanner.scan(tmp_path / FILENAME, nproc=4)
    psf = ParameterScanFile(tmp_path / FILENAME)

    parray0s = psf.parray0s()
    parray1s = psf.parray1s()
    assert psf.parameter_name == PARAMETER_NAME
    parray0 = next(parray0s)
    assert_parray_equality(parray0, pscanner.parray0[0])
    assert_parray_equality(next(parray0s), pscanner.parray0[1])

    assert next(parray1s) != parray0

    assert psf.parameter_name == PARAMETER_NAME

    assert psf.marker_names == set(MARKER_NAMES)
    assert psf.parameter_values == PARAMETER_VALUES
