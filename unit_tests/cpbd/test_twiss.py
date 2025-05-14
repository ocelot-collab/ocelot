import pytest
import numpy as np
from ocelot.cpbd.beam import Twiss
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
    tw = Twiss()
    tw.emit_xn = 0.6e-6
    tw.emit_yn = 0.6e-6
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
    tw = Twiss(emit_x=0.6e-6 / gamma, emit_y=0.6e-6 / gamma)
    tw.E = 14
    tw.emit_x = 0.6e-6 / gamma
    tw.emit_y = 0.6e-6 / gamma
    tw.beta_x = 39.53298538303118
    tw.beta_y = 50.0173089403086
    tw.alpha_x = 1.2073201265148732
    tw.alpha_y = -1.8752915148201614
    _assert_twiss_properties(tw)
