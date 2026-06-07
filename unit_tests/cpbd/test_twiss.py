import numpy as np
from unittest.mock import patch
from ocelot.cpbd.beam import Twiss
from ocelot.cpbd.beam import core as beam_core
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
