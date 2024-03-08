from pathlib import Path

import numpy as np
import pytest

from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.transformations import SecondTM
from ocelot.cpbd.track import ParameterScanner
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd.io import ParameterScanFile
from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.elements import Marker, SBend

PARAMETER_NAME = "pytest-param-name"
PARAMETER_VALUES = [0.3, 0.2]
PARRAY0S = [ParticleArray(10), ParticleArray(20)]
MARKER_NAMES = ["pytest-marker-name-1", "pytest-marker-name-2"]
FILENAME = "test.hdf5"

try:
    import h5py
except ImportError:
    pytest.skip(reason="No h5py installed", allow_module_level=True)

def assert_parray_equality(first, second):
    np.testing.assert_equal(first.E, second.E)
    np.testing.assert_equal(first.s, second.s)
    np.testing.assert_equal(first.q_array, second.q_array)
    np.testing.assert_equal(first.rparticles, second.rparticles)

@pytest.fixture
def pscanner_file(tmp_path):
    f = ParameterScanFile(tmp_path / FILENAME, "w")

    for pval, parray0 in zip(PARAMETER_VALUES, PARRAY0S):
        f.new_run_output(pval, parray0, MARKER_NAMES)
    yield f
    f.close()

def test_ParameterScannerFile_run_names(pscanner_file):
    assert pscanner_file.run_names == ["run-0", "run-1"]

def test_ParameterScannerFile_marker_names(pscanner_file):
    assert pscanner_file.marker_names == {"pytest-marker-name-1", "pytest-marker-name-2"}

def test_ParameterScannerFile_parameter_values(pscanner_file):
    assert pscanner_file.parameter_values == PARAMETER_VALUES

def test_ParameterScannerFile_write_parray0(pscanner_file):
    # Assuming len(parray1) == len(parray0)!
    parray0 = ParticleArray.random(10)
    pscanner_file.write_parray0(0, parray0)

    parray0_read = next(pscanner_file.parray0s())
    assert_parray_equality(parray0, parray0_read)

def test_ParameterScannerFile_write_parray1(pscanner_file):
    # Assuming len(parray1) == len(parray0)!
    parray1 = ParticleArray.random(10)
    pscanner_file.write_parray1(0, parray1)

    parray1_read = next(pscanner_file.parray1s())
    assert_parray_equality(parray1, parray1_read)

def test_ParameterScannerFile_write_marker(pscanner_file):
    marker_name = "pytest-marker-name-1"
    parray_marker = ParticleArray.random(10)
    pscanner_file.write_parray_marker(0, marker_name, parray_marker)

    parray_marker_read = next(pscanner_file.parray_markers(marker_name))
    assert_parray_equality(parray_marker, parray_marker_read)

def test_ParameterScannerFile_set_parameter_name(pscanner_file):
    new_name = "Hello there"
    pscanner_file.parameter_name = new_name
    assert pscanner_file.parameter_name == new_name

def test_ParameterScannerFile_next_run_name(pscanner_file):
    assert pscanner_file.next_run_name() == "run-2"

def test_ParameterScannerFile_filename(pscanner_file):
    assert Path(pscanner_file.filename).name == FILENAME

def test_ParameterScannerFile_init_from_parameter_scanner(tmp_path):

    b1 = SBend(l=0.5, angle=-0.5)
    cell = [Marker(MARKER_NAMES[0]), b1, Marker(MARKER_NAMES[1])]

    magnetic_lattice = MagneticLattice(cell, method={"global": SecondTM})

    with ParameterScanFile(tmp_path / FILENAME, "w") as f:
        pscanner = ParameterScanner(Navigator(magnetic_lattice),
                                    PARAMETER_VALUES,
                                    PARRAY0S,
                                    PARAMETER_NAME,
                                    [Marker(n) for n in MARKER_NAMES]
                                    )

        f.init_from_parameter_scanner(pscanner)

        parray0s = f.parray0s()
        assert_parray_equality(next(parray0s), PARRAY0S[0])
        assert_parray_equality(next(parray0s), PARRAY0S[1])
        assert f.parameter_name == PARAMETER_NAME
        assert f.parameter_values == PARAMETER_VALUES
        assert f.marker_names == set(MARKER_NAMES)
