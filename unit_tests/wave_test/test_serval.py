"""Small-grid checks for the SERVAL generator's public calling convention."""

import numpy as np
import pytest

from ocelot.common.globals import h_eV_s, speed_of_light
from ocelot.optics import wave


@pytest.mark.parametrize("plot_options", [{}, {"showfig": False}, {"showfig": True}])
def test_serval_grid_and_optional_plot(monkeypatch, plot_options):
    from ocelot.gui import dfl_plot

    plotted = []
    monkeypatch.setattr(dfl_plot, "plot_dfl", lambda field, **kwargs: plotted.append(field))
    monkeypatch.setattr(wave, "nthread", 1)
    field = wave.dfl_gen_undulator_serval(
        E_ph=2167, L_w=3.6, shape=(15, 17, 2), dgrid=(400e-6, 300e-6, 1),
        sig_x=38e-6, sig_y=4.68e-6, sig_xp=25e-6, sig_yp=20e-6,
        seed=1, **plot_options,
    )
    assert field.fld.shape == (2, 17, 15)
    assert field.dx == pytest.approx(400e-6 / 15)
    assert field.dy == pytest.approx(300e-6 / 17)
    assert field.dz == pytest.approx(0.5)
    assert field.xlamds == pytest.approx(h_eV_s * speed_of_light / 2167)
    assert np.isfinite(field.fld).all()
    assert np.any(field.fld != 0)
    assert len(plotted) == int(plot_options.get("showfig", False))


def test_zero_emittance_beam_pointing_changes_the_angle_of_each_slice(monkeypatch):
    monkeypatch.setattr(wave, "nthread", 1)
    grid = 300.e-6
    wavelength = h_eV_s * speed_of_light / 2167
    # One transverse Fourier bin corresponds to lambda / window width.
    angle_step = wavelength / grid
    field = wave.dfl_gen_undulator_serval(
        E_ph=2167, shape=(15, 17, 2), dgrid=(grid, grid, 1.),
        jitter_kx=[0., 2 * angle_step], s_support="beam", k_support="beam", seed=42)
    assert np.isfinite(field.fld).all()
    for index, peak_x in enumerate((7, 9)):
        power = np.abs(field.fld[index])**2
        assert power[8, peak_x] > 0
        assert np.count_nonzero(power) == 1


def test_intensity_support_keeps_the_initial_random_field_phase(monkeypatch):
    # Capture the field on the source plane before transforming it to angle.
    source = []
    transform = wave.RadiationField.to_domain

    def capture(self, domains, **kwargs):
        if domains == "kf":
            source.append(self.fld.copy())
        return transform(self, domains, **kwargs)

    monkeypatch.setattr(wave.RadiationField, "to_domain", capture)
    monkeypatch.setattr(wave, "nthread", 1)
    wave.dfl_gen_undulator_serval(shape=(21, 19, 2), sig_x=40e-6, sig_y=40e-6, seed=5)
    rng = np.random.RandomState(5)
    shape = source[0].shape
    noise = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    amplitude = source[0] / noise
    np.testing.assert_allclose(amplitude.imag, 0., atol=1e-15)
    assert amplitude.real.min() >= 0.


@pytest.mark.parametrize("mode", ["amplitude", "intensity"])
@pytest.mark.parametrize("bin_offset", [0., 2., 2.35])
def test_even_serval_grid_matches_point_emitter_on_both_planes(monkeypatch, mode, bin_offset):
    from scipy.special import sici

    source, spectrum = [], []
    transform = wave.RadiationField.to_domain

    def capture(self, domains, **kwargs):
        if domains == "kf":
            source.append(self.fld.copy())
        result = transform(self, domains, method="np")
        if domains == "kf":
            spectrum.append(self.fld.copy())
        return result

    monkeypatch.setattr(wave.RadiationField, "to_domain", capture)
    nx, ny, grid_x, grid_y = 12, 10, 150e-6, 160e-6
    energy, length = 2000., 2.6
    wavelength = h_eV_s * speed_of_light / energy
    pointing = bin_offset * wavelength / grid_x
    field = wave.dfl_gen_undulator_serval(
        E_ph=energy, L_w=length, shape=(nx, ny, 1), dgrid=(grid_x, grid_y, 1.),
        s_support=mode, k_support=mode, jitter_kx=pointing, seed=9)

    rng = np.random.RandomState(9)
    noise = rng.normal(size=(1, ny, nx)) + 1j * rng.normal(size=(1, ny, nx))
    # The ideal electron at x=y=0 radiates symmetrically on the real mesh,
    # whose even aperture has two observation pixels straddling each axis.
    x = np.linspace(-grid_x / 2, grid_x / 2, nx)
    y = np.linspace(-grid_y / 2, grid_y / 2, ny)[:, None]
    near = 1j * (np.pi - 2 * sici(2*np.pi*(x**2 + y**2)/(wavelength * length))[0])
    if mode == "intensity":
        near = abs(near)
    np.testing.assert_allclose(source[0], noise * near, atol=1e-13)

    # In angular space k=0 IS a pixel. The displaced sinc centre must lie
    # at the electron pointing, even when it is between two bins.
    angle_x = (np.arange(nx) - nx//2) * wavelength / grid_x
    angle_y = (np.arange(ny) - ny//2)[:, None] * wavelength / grid_y
    far = np.sinc(length * ((angle_x - pointing)**2 + angle_y**2) / (2*wavelength))
    if mode == "intensity":
        far = abs(far)
    np.testing.assert_allclose(field.fld, spectrum[0] * far / np.sum(far), atol=2e-13)


@pytest.mark.parametrize("nx,ny,mx,my", [(12, 10, 0, 0), (12, 10, 2, -1), (11, 9, -2, 1)])
def test_diffraction_propagates_a_single_fourier_mode_with_its_physical_phase(monkeypatch, nx, ny, mx, my):
    # A tilted plane wave is an eigenfunction of vacuum propagation.
    monkeypatch.setattr(wave, "fftw_avail", False)
    field = wave.RadiationField((1, ny, nx))
    field.dx, field.dy, field.dz, field.xlamds = 3e-4, 2e-4, 1e-3, 4e-5
    plane = np.exp(2j*np.pi*(mx*np.arange(nx)/nx + my*np.arange(ny)[:, None]/ny))
    field.fld[0] = plane
    distance = .03
    field.prop(distance, fine=-1)
    k = 2*np.pi / field.xlamds
    k_perp2 = (2*np.pi*mx/(nx*field.dx))**2 + (2*np.pi*my/(ny*field.dy))**2
    phase = -distance * k_perp2 / (np.sqrt(k**2 - k_perp2) + k)
    np.testing.assert_allclose(field.fld[0], plane*np.exp(1j*phase), atol=1e-12)
