import numpy as np
import pytest


@pytest.fixture
def without_trapz(monkeypatch):
    monkeypatch.delattr(np, "trapz", raising=False)


def test_beam_charge_on_nonuniform_grid(without_trapz):
    from ocelot.common.globals import speed_of_light
    from ocelot.cpbd.beam import BeamArray

    beam = BeamArray(3)
    beam.s = np.array([0.0, 1.0, 3.0]) * 1e-15 * speed_of_light
    beam.I = np.array([2.0, 4.0, 8.0])

    assert beam.charge() == pytest.approx(15e-15, abs=1e-28)


def test_fel_energy_on_nonuniform_grid(without_trapz, monkeypatch):
    from ocelot.common.globals import speed_of_light
    from ocelot.rad.fel import FelParameters

    parameters = FelParameters()
    parameters.s = np.array([0.0, 1.0, 3.0]) * 1e-15 * speed_of_light
    monkeypatch.setattr(parameters, "P", lambda z: np.array([2.0, 4.0, 8.0]))

    assert parameters.E() == pytest.approx(15e-15, abs=1e-28)


@pytest.mark.parametrize("module_name", ["wave", "new_wave"])
@pytest.mark.parametrize("spec_squared", [False, True])
@pytest.mark.parametrize("multiple_events", [False, True])
def test_photon_spectrum_normalization(without_trapz, module_name, spec_squared, multiple_events):
    from importlib import import_module

    module = import_module(f"ocelot.optics.{module_name}")
    energy = np.array([100.0, 101.0, 103.0, 106.0])
    spectrum = np.array([1.0, 2.0, 3.0, 1.0])
    photons = 25.0
    if multiple_events:
        spectrum = np.column_stack([spectrum, 2 * spectrum, np.zeros_like(spectrum)])
        photons = np.array([25.0, 50.0, 0.0])
    if not spec_squared:
        spectrum = spectrum * (1 + 2j)

    result = module.calc_ph_sp_dens(spectrum, energy, photons, spec_squared=spec_squared)
    density = result if spec_squared else abs(result) ** 2
    widths = np.diff(energy)
    if multiple_events:
        widths = widths[:, None]
    integrated = np.sum((density[:-1] + density[1:]) * widths / 2, axis=0)

    assert result.shape == spectrum.shape
    np.testing.assert_allclose(integrated, photons)


@pytest.mark.parametrize("module_name", ["wave", "new_wave"])
def test_sase_mean_pulse_energy_from_envelope(without_trapz, module_name):
    from importlib import import_module
    from ocelot.common.globals import speed_of_light

    module = import_module(f"ocelot.optics.{module_name}")
    positions = np.linspace(0.0, 3e-6, 32)
    envelope = np.full(32, 2e6)
    energies = np.linspace(100.0, 1000.0, 32)
    scale, field, _, _ = module.imitate_1d_sase_like(
        positions, envelope, energies, np.ones(32), n_events=3, seed=7,
    )
    power = abs(field) ** 2
    pulse_energies = np.sum(
        (power[:-1] + power[1:]) * np.diff(scale)[:, None] / (2 * speed_of_light), axis=0,
    )

    assert pulse_energies.mean() == pytest.approx(2e6 * 3e-6 / speed_of_light, rel=1e-12, abs=0)


def test_genesis4_charge_energy_and_spectrum(without_trapz):
    from ocelot.adaptors.genesis4 import Genesis4Output
    from ocelot.common.globals import speed_of_light

    output = Genesis4Output()
    output.tdp = True
    output.s = np.arange(4) * 1e-15 * speed_of_light
    output.h5 = {
        "Beam/current": np.array([[2.0, 4.0, 6.0, 8.0]]),
        "Field/power": np.array([[2.0, 4.0, 6.0, 8.0], [4.0, 8.0, 12.0, 16.0]]),
        "Field/intensity-nearfield": np.ones((2, 4)),
        "Field/phase-nearfield": np.zeros((2, 4)),
        "Global/lambdaref": np.array([1e-9]),
        "Global/s": output.s,
    }

    assert output.beam_charge == pytest.approx(15e-15, abs=1e-28)
    np.testing.assert_allclose(output.rad_energy, [15e-15, 30e-15], atol=0, rtol=1e-12)
    for zi in (None, 1):
        energy, spectrum = output.calc_spec(zi=zi)
        photons = np.sum((spectrum[..., :-1] + spectrum[..., 1:]) * np.diff(energy) / 2, axis=-1)
        expected = output.n_photons if zi is None else output.n_photons[zi]
        np.testing.assert_allclose(photons, expected)


def test_crisp_current_normalization(without_trapz):
    from ocelot.utils.crisp import norm_current

    time = np.array([0.0, 1.0, 3.0]) * 1e-15
    current = norm_current(time, np.array([2.0, 4.0, 8.0]), charge=150e-15)

    np.testing.assert_allclose(current, [20.0, 40.0, 80.0])


def test_math_helpers_use_supported_numpy_dtypes():
    from ocelot.common.math_op import bin_scale, conj_sym

    result = conj_sym([1, 2 + 3j, 3, 4])

    np.testing.assert_array_equal(result, [1, 2 + 3j, 3, 2 - 3j])
    np.testing.assert_array_equal(bin_scale(np.arange(10.0), 2), [1, 3, 5, 7, 9])


def test_response_matrix_set_membership_without_in1d(monkeypatch):
    from ocelot.cpbd.response_matrix import ResponseMatrix, ResponseMatrixJSON

    monkeypatch.delattr(np, "in1d", raising=False)

    response = ResponseMatrixJSON()
    response.cor_names = ["c1", "c2", "c3"]
    response.bpm_names = ["b1", "b2"]
    response.matrix = np.arange(12.0).reshape(4, 3)

    selected = response.extract(["c1", "c3"], ["b2"])
    np.testing.assert_array_equal(selected, [[3, 5], [9, 11]])

    replacement = np.full((2, 2), -1.0)
    response.inject(["c1", "c3"], ["b2"], replacement)
    np.testing.assert_array_equal(response.extract(["c1", "c3"], ["b2"]), replacement)
    response.compare(response)

    dataframe_response = ResponseMatrix()
    dataframe_response.cor_names = response.cor_names
    dataframe_response.bpm_names = response.bpm_names
    dataframe_response.matrix = response.matrix
    dataframe_response.compare(dataframe_response)


def test_fel_ptap_at_zero_returns_nan():
    from ocelot.rad.fel import FelParameters

    parameters = FelParameters()
    parameters.P_sn = np.array([1.0, 2.0])

    assert np.isnan(parameters.Ptap(z=0)).all()


def test_optics_signals_use_complex_dtype():
    from ocelot.optics.utils import Signal, Signal3D

    assert np.issubdtype(Signal().f.dtype, np.complexfloating)
    assert np.issubdtype(Signal3D().f.dtype, np.complexfloating)


def test_cut_lattice_accepts_fractional_cell_count():
    from ocelot.adaptors.genesis import cut_lattice
    from ocelot.cpbd.elements import Drift
    from ocelot.cpbd.magnetic_lattice import MagneticLattice

    lattice = MagneticLattice(tuple(Drift(l=1.0, eid=f"d{i}") for i in range(4)))

    trimmed = cut_lattice(lattice, 1.2, elem_in_cell=1)

    assert [element.id for element in trimmed.sequence] == ["d2", "d3"]
