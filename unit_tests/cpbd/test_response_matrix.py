from functools import partial

import numpy as np
import pytest

from ocelot.common.globals import m_e_GeV, speed_of_light
from ocelot.cpbd.beam import Particle, Twiss
from ocelot.cpbd.elements import Cavity, Drift, Hcor, Matrix, Monitor, Quadrupole, SBend, TWCavity, Vcor
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.optics import UnstableLatticeError
from ocelot.cpbd.orbit_correction import Orbit
from ocelot.cpbd.response_matrix import (
    LinacCavityScanRM, LinacDisperseSimRM, LinacDisperseTmatrixRM, LinacOpticalRM, LinacRmatrixRM,
    LinacSimRM, ResponseMatrix, RingRM,
)
from ocelot.cpbd.transformations import ExactDriftTM, SecondTM, TransferMap


@pytest.fixture
def orbit():
    lattice = MagneticLattice((
        Cavity(l=1, v=0.57, phi=0, freq=0),
        Hcor(angle=0.0005, eid="H"), Drift(l=2), Monitor(eid="B1"),
        Vcor(angle=-0.0003, eid="V"), Drift(l=3), Monitor(eid="B2"),
    ), method={"global": SecondTM})
    return Orbit(lattice, rm_method=LinacSimRM, disp_rm_method=LinacDisperseSimRM)


def test_tracking_response_at_nonzero_settings(orbit):
    seed = Twiss(E=0.13)
    orbit.response_matrix.calculate(seed)
    orbit.disp_response_matrix.calculate(seed)

    # Thin correctors in a drift give +L on momentum and -L/beta per
    # unit *local* delta. The on-crest cavity scales entrance delta by p0/p1.
    orm = np.array([[2, 0], [5, 0], [0, 0], [0, 3]])
    p0, p1 = np.sqrt(np.array([seed.E, 0.70])**2 - m_e_GeV**2)
    drm = -orm * (0.70 / p1) * (p0 / p1)
    np.testing.assert_allclose(orbit.response_matrix.matrix, orm, atol=1e-12)
    np.testing.assert_allclose(orbit.disp_response_matrix.matrix, drm, atol=1e-12)

    angles = [orbit.hcors[0].angle, orbit.vcors[0].angle]
    np.testing.assert_array_equal(angles, [0.0005, -0.0003])
    np.testing.assert_allclose(orbit.get_orbit(), orm @ angles, atol=1e-14)
    np.testing.assert_allclose(orbit.get_dispersion(), drm @ angles, atol=1e-14)


@pytest.mark.parametrize("matrix_name", ["response_matrix", "disp_response_matrix"])
def test_failed_probe_restores_correctors(orbit, matrix_name, monkeypatch):
    matrix = getattr(orbit, matrix_name)
    read_orbit = matrix.method.read_virtual_orbit
    hcor, vcor = orbit.hcors[0], orbit.vcors[0]
    angles = (hcor.angle, vcor.angle)

    def fail_on_vertical_probe(*args, **kwargs):
        result = read_orbit(*args, **kwargs)
        if vcor.angle != angles[1]:
            raise RuntimeError("failed BPM read")
        return result

    monkeypatch.setattr(matrix.method, "read_virtual_orbit", fail_on_vertical_probe)
    with pytest.raises(RuntimeError, match="failed BPM read"):
        matrix.calculate(Twiss(E=0.13))

    assert (hcor.angle, vcor.angle) == angles
    np.testing.assert_allclose(orbit.bpms[-1].y, 3 * vcor.angle)
    np.testing.assert_allclose(vcor.tms[0].get_params(0.70).B[3, 0], vcor.angle)


@pytest.mark.parametrize("calculator", [LinacRmatrixRM, LinacOpticalRM, LinacSimRM])
def test_short_line_samples_device_centres(calculator):
    h = Hcor(l=0.08)
    before, after = Monitor(eid="B"), Monitor(l=0.08, eid="B")
    lat = MagneticLattice((before, h, Drift(l=0.03), after))
    # Caller ordering and occurrence identity matter even when IDs coincide.
    method = calculator(lat, [h], [], [after, before])
    np.testing.assert_allclose(method.calculate(Twiss(E=0.13)), [[0.11], [0.], [0.], [0.]], atol=1e-14)
    if calculator is LinacOpticalRM:
        assert after.beta_x > 0.


@pytest.mark.parametrize("calculator", [LinacRmatrixRM, LinacOpticalRM])
def test_coupled_orm_matches_tracking(calculator):
    h, v, bpm = Hcor(l=0.1), Vcor(), Monitor(l=0.12)
    h.tilt = 0.15
    lat = MagneticLattice((h, Drift(l=0.4), v, Quadrupole(l=0.3, k1=1.5, tilt=0.6),
                           Drift(l=1.2), bpm))
    orbit = Orbit(lat, rm_method=LinacSimRM)
    seed = Twiss(E=0.13)
    orbit.response_matrix.calculate(seed)
    expected = orbit.response_matrix.matrix
    assert min(abs(expected[1, 0]), abs(expected[0, 1])) > 0.01
    np.testing.assert_allclose(calculator(lat, [h], [v], [bpm]).calculate(seed), expected, atol=1e-13)


@pytest.mark.parametrize("drift_tm", [SecondTM, TransferMap, ExactDriftTM])
def test_direct_drm_uses_active_map_and_entrance_normalization(drift_tm):
    h, v, bpm = Hcor(), Vcor(), Monitor()
    lat = MagneticLattice((Cavity(l=1, v=0.57), h, v, Drift(l=2), bpm), method={Drift: drift_tm})
    seed = Twiss(E=0.13)
    matrix = LinacDisperseTmatrixRM(lat, [h], [v], [bpm]).calculate(seed)
    p0, p1 = np.sqrt(np.array([seed.E, 0.70])**2 - m_e_GeV**2)
    slope = 0. if drift_tm is TransferMap else -2. * (0.70 / p1) * (p0 / p1)
    np.testing.assert_allclose(matrix, slope * np.eye(2), atol=2e-9)


@pytest.mark.parametrize("cavity", [Cavity, TWCavity])
def test_direct_drm_matches_mixed_tracking_difference(cavity):
    h, v = Hcor(l=0.1, angle=0.0007), Vcor(angle=-0.0003)
    h.tilt = 0.2
    before, middle, after = Monitor(), Monitor(), Monitor()
    lat = MagneticLattice((
        before, SBend(l=0.35, angle=0.09), h, Drift(l=0.5), middle,
        cavity(l=0.8, v=0.08, freq=1.3e9, phi=28), v,
        Quadrupole(l=0.25, k1=1.4, k2=2., tilt=0.42),
        SBend(l=0.4, angle=-0.13), Drift(l=1.5), after,
    ), method={"global": SecondTM})
    seed = Twiss(E=0.13)
    original = (h.angle, v.angle)
    method = LinacDisperseTmatrixRM(lat, [h], [v], [before, middle, after])
    result = method.calculate(seed)
    assert (h.angle, v.angle) == original

    def track_bpms(delta):
        particle = Particle(E=seed.E, p=delta)
        readings = []
        for element in lat.sequence:
            for tm in element.tms:
                tm.apply(particle)
            if isinstance(element, Monitor):
                readings.append((particle.x, particle.y))
        return np.array(readings).T.ravel()

    dtheta, dp = 1e-5, 1e-5
    for j, cor in enumerate((h, v)):
        cor.angle = original[j] + dtheta
        plus = track_bpms(dp) - track_bpms(-dp)
        cor.angle = original[j] - dtheta
        minus = track_bpms(dp) - track_bpms(-dp)
        cor.angle = original[j]
        expected = (plus - minus) / (4 * dtheta * dp)
        np.testing.assert_allclose(result[:, j], expected, rtol=2e-6, atol=2e-8)
    assert result[0, 0] == result[1, 1] == 0.
    assert abs(result[5, 0]) > 0.01


@pytest.mark.parametrize("cavity", [Cavity, TWCavity])
def test_rf_curvature_at_the_actual_arrival_phase(cavity):
    # tau = a*theta + b*delta before the cavity. Read outgoing delta as x.
    # The mixed derivative isolates the curvature of the RF sinusoid.
    h, bpm = Hcor(angle=0.013), Monitor()
    a, b = 0.3, -0.25
    before, after = Matrix(), Matrix()
    before.r = np.eye(6)
    before.r[4, 1], before.r[4, 5] = a, b
    after.r = np.eye(6)
    after.r[0] = after.r[5]
    rf = cavity(l=0.9, v=0.04, freq=1.3e9, phi=32)
    seed = Twiss(E=0.13)
    lat = MagneticLattice((h, before, rf, after, bpm))

    phi = np.deg2rad(rf.phi)
    k = 2 * np.pi * rf.freq / speed_of_light
    pc_out = np.sqrt((seed.E + rf.v * np.cos(phi))**2 - m_e_GeV**2)
    expected = -a * b * rf.v * k**2 * np.cos(phi - k * a * h.angle) / pc_out
    result = LinacDisperseTmatrixRM(lat, [h], [], [bpm]).calculate(seed)
    np.testing.assert_allclose(result[:, 0], [expected, 0.], atol=1e-13)


@pytest.mark.parametrize("cavity", [Cavity, TWCavity])
@pytest.mark.parametrize("step", [0.001, -0.002])
def test_cavity_scan_fixed_magnet_response(cavity, step):
    rf = cavity(l=0.8, v=0.035, freq=1.3e9, phi=33.)
    h, v = Hcor(angle=0.0007, eid="H"), Vcor(angle=-0.0004, eid="V")
    upstream, downstream = Monitor(eid="B1"), Monitor(l=0.2, eid="B2")
    lat = MagneticLattice((upstream, rf, h, v, Drift(l=1.4), downstream), method={"global": SecondTM})
    energy, normalization = 0.13, 0.7
    scan = partial(LinacCavityScanRM, cavity_weights={rf: 1.}, voltage_step=step,
                   normalization_energy=normalization)
    orbit = Orbit(lat, disp_rm_method=scan)
    rm = orbit.disp_response_matrix
    maps = rf.tms
    rm.calculate(Twiss(E=energy))

    # The scan injects dV*cos(phi)/pc at the correctors. In a second-order
    # drift the slope per unit local delta is -L/beta for either kick.
    cosine = np.cos(np.deg2rad(rf.phi))
    e1 = energy + rf.v * cosine
    pc = np.sqrt(e1**2 - m_e_GeV**2)
    length = 1.4 + downstream.l / 2.
    slope = -length * (e1 / pc) * normalization * cosine / pc
    expected = slope * np.array([[0, 0], [1, 0], [0, 0], [0, 1]])
    np.testing.assert_allclose(rm.get_matrix(), expected, atol=3e-12)
    np.testing.assert_array_equal(rm.matrix[[0, 2]], 0.)
    np.testing.assert_allclose(orbit.get_dispersion(), expected @ [h.angle, v.angle], atol=1e-14)
    np.testing.assert_allclose(orbit.get_orbit(), [0, length * h.angle, 0, length * v.angle])
    np.testing.assert_allclose([upstream.E, downstream.E], [energy, e1])
    assert (rf.v, h.angle, v.angle) == (0.035, 0.0007, -0.0004)
    assert rf.tms is maps


def test_station_scan_energy_and_arrival_phase():
    first = Cavity(l=1.3, v=0.0009, phi=12., eid="RF")
    second = TWCavity(l=0.7, v=0.0008, phi=32., freq=1.3e9, eid="RF")
    monitor = Monitor()
    diagnostic = Matrix()
    diagnostic.r = np.eye(6)
    diagnostic.r[0] = diagnostic.r[5]  # Measure outgoing particle delta as X.
    lat = MagneticLattice((first, second, diagnostic, monitor))
    e0, step, normalization = 0.0015, 0.00011, 0.004
    scan = LinacCavityScanRM(lat, [], [], [monitor], cavity_weights={second: 0.75, first: 0.25},
                             voltage_step=step, normalization_energy=normalization)
    monitor.Dx = 99.
    dx, dy = scan.read_virtual_dispersion(e0, write2bpms=False)
    assert monitor.Dx == 99.

    def pc(energy):
        return np.sqrt(energy**2 - m_e_GeV**2)

    # Absolute first-cavity transit time and the energy gained at the second
    # cavity's actual arrival phase, independently of Ocelot's coordinate change.
    phi1, phi2 = np.deg2rad([first.phi, second.phi])
    e1 = e0 + first.v * np.cos(phi1)
    e1_scan = e0 + (first.v + 0.25 * step) * np.cos(phi1)
    tau = first.l * ((e0 + e1_scan) / (pc(e0) + pc(e1_scan)) - (e0 + e1) / (pc(e0) + pc(e1)))
    e2 = e1 + second.v * np.cos(phi2)
    phase = phi2 - 2 * np.pi * second.freq * tau / speed_of_light
    e2_scan = e1_scan + (second.v + 0.75 * step) * np.cos(phase)
    expected = (e2_scan - e2) / pc(e2) / (step / normalization)
    np.testing.assert_allclose(dx, expected, rtol=2e-12)
    np.testing.assert_array_equal(dy, [0.])


def test_cavity_scan_corrector_before_station_matches_finite_kicks():
    h, v = Hcor(angle=0.002), Vcor(angle=-0.001)
    first = Cavity(l=0.7, v=0.035, phi=15., freq=1.3e9, vx_up=0.0008)
    second = TWCavity(l=0.5, v=0.04, phi=28., freq=1.3e9)
    b0, b1, b2 = Monitor(), Monitor(), Monitor()
    lat = MagneticLattice((h, v, Drift(l=1.2), b0, first, SBend(l=0.4, angle=0.09),
                           b1, second, Quadrupole(l=0.3, k1=1.4, tilt=0.3), Drift(l=1.), b2),
                          method={"global": SecondTM})
    scan = LinacCavityScanRM(lat, [h], [v], [b2, b0, b1], cavity_weights={first: 0.4, second: 0.6},
                             voltage_step=0.002, normalization_energy=0.25)
    seed = Twiss(E=0.13)
    matrix = scan.calculate(seed)
    np.testing.assert_array_equal(matrix[[1, 4]], 0.)
    assert abs(matrix[0, 1]) > 0.01  # rolled quadrupole couples Y kicks into X
    for i, corrector in enumerate((h, v)):
        angle, step = corrector.angle, 1e-5
        corrector.angle = angle + step
        plus = np.concatenate(scan.read_virtual_dispersion(seed.E, write2bpms=False))
        corrector.angle = angle - step
        minus = np.concatenate(scan.read_virtual_dispersion(seed.E, write2bpms=False))
        corrector.angle = angle
        np.testing.assert_allclose(matrix[:, i], (plus - minus) / (2 * step), atol=1e-8)


def test_invalid_cavity_scan_leaves_lattice_and_bpm_readings_intact():
    rf, bpm = Cavity(l=0.5, v=0.02), Monitor()
    lat = MagneticLattice((rf, bpm))
    options = dict(cavity_weights={rf: 1.}, voltage_step=0.001, normalization_energy=0.7)
    for bad in (dict(voltage_step=0.), dict(voltage_step=np.nan), dict(normalization_energy=-1.),
                dict(cavity_weights={}), dict(cavity_weights={rf: np.inf}), dict(cavity_weights={rf: 0.})):
        with pytest.raises(ValueError):
            LinacCavityScanRM(lat, [], [], [bpm], **{**options, **bad})
    with pytest.raises(ValueError, match="not present"):
        LinacCavityScanRM(lat, [], [], [bpm], **{**options, "cavity_weights": {Cavity(): 1.}})
    with pytest.raises(ValueError, match="appears 2 times"):
        LinacCavityScanRM(MagneticLattice((rf, rf, bpm)), [], [], [bpm], **options)
    with pytest.raises(TypeError, match="Cavity"):
        LinacCavityScanRM(lat, [], [], [bpm], **{**options, "cavity_weights": {bpm: 1.}})
    scan = LinacCavityScanRM(lat, [], [], [bpm], **{**options, "voltage_step": -0.2})
    with pytest.raises(ValueError, match="entrance energy"):
        scan.calculate()
    bpm.Dx, maps = 99., rf.tms
    with pytest.raises(ValueError, match="rest energy"):
        scan.calculate(Twiss(E=0.13))
    assert rf.v == 0.02 and rf.tms is maps and bpm.Dx == 99.


def test_coupled_ring_response_closes_after_one_turn():
    bpm, h, v = Monitor(), Hcor(l=0.1), Vcor()
    lat = MagneticLattice((bpm, Drift(l=0.65), h, Drift(l=0.2), Quadrupole(l=0.25, k1=1.1),
                           v, Drift(l=1.1), Quadrupole(l=0.25, k1=-1.05, tilt=0.005), Drift(l=0.55)))
    response = RingRM(lat, [h], [v], [bpm]).calculate(Twiss(E=0.13))

    def one_turn(launch):
        p = Particle(x=launch[0], px=launch[1], y=launch[2], py=launch[3], E=0.13)
        for elem in lat.sequence:
            for tm in elem.tms:
                tm.apply(p)
        return np.array([p.x, p.px, p.y, p.py])

    R = np.column_stack([one_turn(basis) for basis in np.eye(4)])
    for j, cor in enumerate((h, v)):
        cor.angle = 1e-5
        launch = np.linalg.solve(np.eye(4) - R, one_turn(np.zeros(4)))
        np.testing.assert_allclose(one_turn(launch), launch, atol=1e-15)
        np.testing.assert_allclose(response[:, j], launch[[0, 2]] / cor.angle, rtol=1e-10)
        cor.angle = 0.
    assert abs(response[1, 0]) > 0.1


def test_drift_has_no_periodic_ring_response():
    h, bpm = Hcor(), Monitor()
    lat = MagneticLattice((h, Drift(l=1), bpm))
    with pytest.raises(UnstableLatticeError):
        RingRM(lat, [h], [], [bpm]).calculate()


@pytest.mark.parametrize("device", [Hcor, Monitor])
@pytest.mark.parametrize("calculator", [LinacRmatrixRM, LinacSimRM, LinacDisperseTmatrixRM])
def test_ambiguous_response_device_is_rejected(device, calculator):
    repeat = device()
    h = repeat if device is Hcor else Hcor()
    bpm = repeat if device is Monitor else Monitor()
    lat = MagneticLattice((h, Drift(l=1), bpm, repeat))
    orbit = Orbit(lat, rm_method=calculator)
    with pytest.raises(ValueError, match="appears 2 times"):
        orbit.response_matrix.calculate(Twiss(E=0.13))


def test_fitted_matrix_matches_extract_and_save(tmp_path):
    import pandas as pd

    matrix = ResponseMatrix()
    matrix.cor_names, matrix.bpm_names = ["H", "V"], ["B"]
    kicks = np.array([[0, 0], [1, 0], [0, 1], [1, 1]])
    expected = np.array([[2., 0.4], [-0.3, 1.5]])
    data = np.column_stack((kicks, kicks @ expected.T + [2., 1.]))
    fitted = matrix.retrieve_from_scan(pd.DataFrame(data, columns=["H", "V", "B.X", "B.Y"]))
    np.testing.assert_allclose(fitted, expected)
    np.testing.assert_allclose(matrix.get_matrix(), expected)
    np.testing.assert_allclose(matrix.extract(["V"], ["B"]), expected[:, 1:])
    filename = tmp_path / "response.p"
    matrix.dump(filename)
    loaded = ResponseMatrix()
    loaded.load(filename)
    np.testing.assert_allclose(loaded.get_matrix(), expected)


@pytest.mark.parametrize("calculator", [LinacRmatrixRM, RingRM])
def test_clean_matrix_preserves_ring_response_and_updates_array(calculator):
    upstream, downstream = Monitor(eid="B1"), Monitor(eid="B2")
    h = Hcor(eid="H")
    lat = MagneticLattice((upstream, Drift(l=1), h, Drift(l=1), downstream))
    orbit = Orbit(lat, rm_method=calculator)
    matrix = orbit.response_matrix
    matrix.cor_names, matrix.bpm_names = ["H"], ["B1", "B2"]
    matrix.df = matrix.data2df(np.ones((4, 1)), matrix.bpm_names, matrix.cor_names)
    matrix.clean_rm(coupling=True)
    expected = [0., 1., 0., 1.] if calculator is LinacRmatrixRM else [1.] * 4
    np.testing.assert_array_equal(matrix.get_matrix()[:, 0], expected)
    matrix.clean_rm(coupling=False)
    np.testing.assert_array_equal(matrix.get_matrix()[:, 0], [*expected[:2], 0., 0.])
