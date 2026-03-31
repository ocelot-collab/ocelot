from pathlib import Path

import numpy as np

from ocelot.adaptors.elegant_lattice_converter import ElegantLatticeConverter
from ocelot.cpbd.elements import Cavity, Matrix, Undulator
from ocelot.cpbd.magnetic_lattice import MagneticLattice


def test_elegant2ocelot_uses_public_wrapper_api_for_parameter_assignment(tmp_path):
    lattice_file = tmp_path / "wrapper_input.lte"
    lattice_file.write_text(
        "M1: EMATRIX, L=0.4, R11=1.2, R56=0.03\n"
        "C1: RFCA, L=0.5, VOLT=20000000, FREQ=1300000000, PHASE=60\n"
        "U1: WIGGLER, L=2.4, K=1.1, POLES=12\n"
        "CELL: LINE=(M1, C1, U1)\n"
        "USE,CELL\n"
        "RETURN\n"
    )

    matrix, cavity, undulator = ElegantLatticeConverter().elegant2ocelot(str(lattice_file))

    assert isinstance(matrix, Matrix)
    assert np.isclose(matrix.l, 0.4)
    assert np.isclose(matrix.r[0, 0], 1.2)
    assert np.isclose(matrix.r[4, 5], 0.03)

    assert isinstance(cavity, Cavity)
    assert np.isclose(cavity.l, 0.5)
    assert np.isclose(cavity.v, 0.02)
    assert np.isclose(cavity.freq, 1.3e9)
    assert np.isclose(cavity.phi, 30.0)

    assert isinstance(undulator, Undulator)
    assert np.isclose(undulator.l, 2.4)
    assert np.isclose(undulator.Kx, 1.1)
    assert np.isclose(undulator.nperiods, 6.0)
    assert np.isclose(undulator.lperiod, 0.4)


def test_ocelot2elegant_uses_public_wrapper_api_for_parameter_export(tmp_path):
    matrix = Matrix(l=0.4, eid="M1")
    matrix.r[0, 0] = 1.2
    matrix.r[4, 5] = 0.03
    cavity = Cavity(l=0.5, v=0.02, freq=1.3e9, phi=30.0, eid="C1")
    lattice = MagneticLattice((matrix, cavity))
    output_file = tmp_path / "wrapper_output.lte"

    ElegantLatticeConverter().ocelot2elegant(lattice, file_name=str(output_file))
    output = Path(output_file).read_text()

    assert "M1: EMATRIX,L=0.4" in output
    assert ",R11=1.2" in output
    assert ",R56=0.03" in output
    assert "C1: RFCA,L=0.5" in output
    assert ",VOLT=20000000.0" in output
    assert ",FREQ=1300000000.0" in output
    assert ",PHASE=60.0" in output
