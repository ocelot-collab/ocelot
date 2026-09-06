import numpy as np
import pytest
from unittest.mock import patch
from ocelot.cpbd.beam import Twiss
from ocelot.cpbd.beam import core as beam_core
from ocelot.cpbd.elements import Drift, Quadrupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.optics import UnstableLatticeError, periodic_twiss, twiss
from ocelot.common.globals import m_e_GeV

gamma = 14 / m_e_GeV

def _assert_twiss_properties(tw):
    assert np.isclose(tw.emit_xn, 0.6e-6, rtol=1e-10)
    assert np.isclose(tw.emit_yn, 0.6e-6, rtol=1e-10)
    assert np.isclose(tw.beta_x, 39.53298538303118, rtol=1e-10)
    assert np.isclose(tw.beta_y, 50.0173089403086, rtol=1e-10)
    assert np.isclose(tw.alpha_x, 1.2073201265148732, rtol=1e-10)
    assert np.isclose(tw.alpha_y, -1.8752915148201614, rtol=1e-10)
    assert np.isclose(tw.gamma_x, 0.06216636219289371, rtol=1e-10)
    assert np.isclose(tw.gamma_y, 0.09030310429029308, rtol=1e-10)
    assert np.isclose(tw.E, 14, rtol=1e-10)

def test_tws1():
    tw = Twiss(
        emit_x=0.6e-6 / gamma,
        emit_y=0.6e-6 / gamma,
        beta_x=39.53298538303118,
        beta_y=50.0173089403086,
        alpha_x=1.2073201265148732,
        alpha_y=-1.8752915148201614,
        E=14,
    )
    _assert_twiss_properties(tw)

def test_tws2():
    tw = Twiss(
        emit_xn=0.6e-6,
        emit_yn=0.6e-6,
        beta_x=39.53298538303118,
        beta_y=50.0173089403086,
        alpha_x=1.2073201265148732,
        alpha_y=-1.8752915148201614,
        E=14,
    )
    _assert_twiss_properties(tw)

def test_tws3():
    with patch.object(beam_core.warnings, "warn") as warn:
        tw = Twiss()
        tw.emit_xn = 0.6e-6
        tw.emit_yn = 0.6e-6

    assert warn.call_count == 2
    messages = [call.args[0] for call in warn.call_args_list]
    assert any("Twiss.emit_xn was set while E is 0.0 GeV" in message for message in messages)
    assert any("Twiss.emit_yn was set while E is 0.0 GeV" in message for message in messages)

    tw.beta_x = 39.53298538303118
    tw.beta_y = 50.0173089403086
    tw.alpha_x = 1.2073201265148732
    tw.alpha_y = -1.8752915148201614
    tw.E = 14
    _assert_twiss_properties(tw)

def test_tws4():
    tw = Twiss()
    tw.emit_x = 0.6e-6 / gamma
    tw.emit_y = 0.6e-6 / gamma
    tw.beta_x = 39.53298538303118
    tw.beta_y = 50.0173089403086
    tw.alpha_x = 1.2073201265148732
    tw.alpha_y = -1.8752915148201614
    tw.E = 14
    _assert_twiss_properties(tw)

def test_tws5():
    tw = Twiss()
    tw.E = 14
    tw.emit_x = 0.6e-6 / gamma
    tw.emit_y = 0.6e-6 / gamma
    tw.beta_x = 39.53298538303118
    tw.beta_y = 50.0173089403086
    tw.alpha_x = 1.2073201265148732
    tw.alpha_y = -1.8752915148201614
    _assert_twiss_properties(tw)

def test_tws6():
    with patch.object(beam_core.warnings, "warn") as warn:
        tw = Twiss(emit_x=0.6e-6 / gamma, emit_y=0.6e-6 / gamma)

    warn.assert_not_called()
    assert np.isclose(tw.emit_x, 0.6e-6 / gamma, rtol=1e-10)
    assert np.isclose(tw.emit_y, 0.6e-6 / gamma, rtol=1e-10)
    assert tw.emit_xn == 0.0
    assert tw.emit_yn == 0.0

    tw.E = 14
    tw.beta_x = 39.53298538303118
    tw.beta_y = 50.0173089403086
    tw.alpha_x = 1.2073201265148732
    tw.alpha_y = -1.8752915148201614
    _assert_twiss_properties(tw)


def test_tws_emit_x_without_energy_is_retained_and_normalized_after_energy_is_set():
    emit_x = 0.6e-6 / gamma

    with patch.object(beam_core.warnings, "warn") as warn:
        tw = Twiss(
            emit_x=emit_x,
            beta_x=10,
            beta_y=12,
            alpha_x=0.2,
            alpha_y=-0.3,
        )

    warn.assert_not_called()
    assert tw.emit_x == emit_x
    assert tw.emit_xn == 0.0

    tw.E = 14

    assert np.isclose(tw.emit_x, emit_x, rtol=1e-10)
    assert np.isclose(tw.emit_xn, 0.6e-6, rtol=1e-10)


def test_tws_emit_xn_without_energy_warns_and_is_applied_after_energy_is_set():
    with patch.object(beam_core.warnings, "warn") as warn:
        tw = Twiss(
            emit_xn=0.6e-6,
            beta_x=10,
            beta_y=12,
            alpha_x=0.2,
            alpha_y=-0.3,
        )

    warn.assert_called_once()
    assert "Twiss.emit_xn was set while E is 0.0 GeV" in warn.call_args.args[0]
    assert tw.emit_xn == 0.6e-6
    assert tw.emit_x == 0.0

    tw.E = 14

    assert np.isclose(tw.emit_xn, 0.6e-6, rtol=1e-10)
    assert np.isclose(tw.emit_x, 0.6e-6 / gamma, rtol=1e-10)


def _optics_seed():
    return Twiss(beta_x=10.0, beta_y=12.0, alpha_x=0.2, alpha_y=-0.3, E=1.0)


def _stable_periodic_lattice():
    d = Drift(l=0.5)
    qf = Quadrupole(l=0.2, k1=0.3)
    qdh = Quadrupole(l=0.1, k1=-0.3)
    return MagneticLattice((qdh, d, d, qf, d, d, qdh))


def test_twiss_requires_explicit_initial_parameters():
    lattice = MagneticLattice((Drift(l=1.0),))

    with pytest.raises(TypeError, match="required positional argument: 'tws0'"):
        twiss(lattice)


def test_explicit_twiss_api_is_available_from_root_facade():
    import ocelot as ocl

    assert ocl.twiss is twiss
    assert ocl.periodic_twiss is periodic_twiss
    assert ocl.UnstableLatticeError is UnstableLatticeError


def test_twiss_does_not_treat_blank_twiss_as_periodic_request():
    lattice = _stable_periodic_lattice()

    with pytest.raises(ValueError, match=r"use periodic_twiss\(lattice"):
        twiss(lattice, Twiss())


def test_periodic_twiss_returns_propagated_periodic_optics():
    lattice = _stable_periodic_lattice()

    result = periodic_twiss(lattice, Twiss(E=0.005))

    assert len(result) == len(lattice.sequence) + 1
    assert result[0].beta_x > 0
    assert result[0].beta_y > 0
    np.testing.assert_allclose(result[-1].beta_x, result[0].beta_x)
    np.testing.assert_allclose(result[-1].beta_y, result[0].beta_y)
    np.testing.assert_allclose(result[-1].alpha_x, result[0].alpha_x, atol=1e-14)
    np.testing.assert_allclose(result[-1].alpha_y, result[0].alpha_y, atol=1e-14)


def test_unstable_periodic_twiss_raises_without_logging(caplog):
    lattice = MagneticLattice((Drift(l=1.0),))

    with caplog.at_level("WARNING", logger="ocelot.cpbd.optics"):
        with pytest.raises(UnstableLatticeError, match=r"x: \|Tr\(R\)/2\|=1") as exc_info:
            periodic_twiss(lattice)

    assert exc_info.value.cos_mu_x == pytest.approx(1.0)
    assert exc_info.value.cos_mu_y == pytest.approx(1.0)
    assert not caplog.records


def test_attach2elem_raises_before_attaching_when_an_element_instance_is_repeated():
    drift = Drift(l=1.0, eid="D")
    quad = Quadrupole(l=0.2, k1=1.0, eid="Q")
    lattice = MagneticLattice((drift, quad, drift))

    with pytest.raises(ValueError, match=r"Element 'D'.*indices \[0, 2\]"):
        twiss(lattice, _optics_seed(), attach2elem=True)

    assert not hasattr(drift, "tws")
    assert not hasattr(quad, "tws")


def test_attach2elem_iterable_attaches_only_selected_unique_elements():
    shared_drift = Drift(l=1.0, eid="D")
    quad = Quadrupole(l=0.2, k1=1.0, eid="Q")
    lattice = MagneticLattice((shared_drift, quad, shared_drift))

    result = twiss(lattice, _optics_seed(), attach2elem=[quad])

    assert quad.tws is result[2]
    assert np.isclose(quad.tws.s, 1.2)
    assert not hasattr(shared_drift, "tws")


def test_attach2elem_iterable_rejects_a_selected_repeated_element():
    shared_drift = Drift(l=1.0, eid="D")
    quad = Quadrupole(l=0.2, k1=1.0, eid="Q")
    lattice = MagneticLattice((shared_drift, quad, shared_drift))

    with pytest.raises(ValueError, match=r"Element 'D'.*indices \[0, 2\]"):
        twiss(lattice, _optics_seed(), attach2elem=[shared_drift])


def test_attach2elem_distinguishes_instances_even_when_ids_match():
    first = Drift(l=1.0, eid="D")
    second = Drift(l=2.0, eid="D")
    lattice = MagneticLattice((first, second))

    twiss(lattice, _optics_seed(), attach2elem=True)

    assert np.isclose(first.tws.s, 1.0)
    assert np.isclose(second.tws.s, 3.0)
    assert first.tws is not second.tws


def test_attach2elem_rejects_sampled_twiss_output():
    drift = Drift(l=1.0, eid="D")
    lattice = MagneticLattice((drift,))

    with pytest.raises(ValueError, match="requires nPoints=None"):
        twiss(lattice, _optics_seed(), nPoints=10, attach2elem=True)
