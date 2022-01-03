import warnings

import pandas as pd
import pytest

from ocelot.adaptors.tfs import optics_from_tfs, MADXLatticeConverter
from ocelot.common.globals import m_e_GeV

MISSING_TFS_DEPENDENCY = False
REASON = "Missing optional dependency TFS-Pandas"
try:
    import tfs
except ImportError:
    MISSING_TFS_DEPENDENCY = True


@pytest.fixture
def optics_df():
    # The set of initial optics keys needed for the conversion

    initial = {"ALFX": 0.6,
               "ALFY": 0.7,
               "BETX": 0.8,
               "BETY": 0.9,
               "DX": 1.0,
               "DY": 1.1,
               "DPX": 1.2,
               "DPY": 1.3,
               "X": 1.4,
               "Y": 1.5,
               "PX": 1.6,
               "PY": 1.7}

    out = tfs.TfsDataFrame([initial])

    # Initial energy, emittances and rel gamma stored in .headers attribute.
    energy =  14
    relgamma = energy / m_e_GeV
    # Technically shouldn't set the attr here but it's ok for testing
    out.headers = {"ENERGY": energy,
                   "EX": 1e-10,
                   "EY": 1e-10,
                   "GAMMA": relgamma}


    return out

@pytest.mark.skipif(MISSING_TFS_DEPENDENCY, reason=REASON)
def test_optics_from_tfs(optics_df):
    twiss0 = optics_from_tfs(optics_df)

    # Check all the DataFrame keys are equal in the resulting Twiss instance
    twiss_df_keys_to_twiss_instance_keys = {"ALFX": "alpha_x",
                                            "ALFY": "alpha_y",
                                            "BETX": "beta_x",
                                            "BETY": "beta_y",
                                            "DX": "Dx",
                                            "DY": "Dy",
                                            "DPX": "Dxp",
                                            "DPY": "Dyp",
                                            "X": "x",
                                            "Y": "y",
                                            "PX": "xp",
                                            "PY": "yp",
                                            }

    for twiss_df_key, twiss_key in twiss_df_keys_to_twiss_instance_keys.items():
        twiss0_value = getattr(twiss0, twiss_key)
        twiss_df_value = optics_df.iloc[0][twiss_df_key]

        assert twiss0_value == twiss_df_value


    # Check the info stored in the header (emittances and beam energy) are also correct.
    relgamma = twiss0.E / m_e_GeV
    emit_x = optics_df.headers["EX"]
    emit_y = optics_df.headers["EY"]
    emit_xn = emit_x * relgamma
    emit_yn = emit_y * relgamma

    assert emit_x == twiss0.emit_x
    assert emit_x == twiss0.emit_y
    assert emit_xn == twiss0.emit_xn
    assert emit_yn == twiss0.emit_yn
    assert optics_df.headers["ENERGY"] == twiss0.E


@pytest.fixture
def converter():
    return MADXLatticeConverter()

@pytest.fixture
def ns():
    """trivial object to hold a namespace, comes with a length and a name"""
    class NS: pass
    namespace = NS()
    namespace.NAME = "hello"
    namespace.L = 3.0
    return namespace

@pytest.mark.skipif(MISSING_TFS_DEPENDENCY, reason=REASON)
def test_MADXLatticeConverter_make_marker(converter, ns):
    ns.NAME = "Hello"
    marker = converter.make_marker(ns)

    assert marker.id == ns.NAME

@pytest.mark.skipif(MISSING_TFS_DEPENDENCY, reason=REASON)
def test_MADXLatticeConverter_make_monitor(converter, ns):
    monitor = converter.make_monitor(ns)
    assert monitor.id == ns.NAME
    assert monitor.l == ns.L

@pytest.mark.skipif(MISSING_TFS_DEPENDENCY, reason=REASON)
def test_MADXLatticeConverter_make_drift(converter, ns):
    drift = converter.make_drift(ns)
    assert drift.id == ns.NAME
    assert drift.l == ns.L

@pytest.mark.skipif(MISSING_TFS_DEPENDENCY, reason=REASON)
def test_MADXLatticeConverter_make_sbend(converter, ns):
    ns.ANGLE = 0.5
    ns.TILT = 0.25
    ns.E1 = +0.1
    ns.E2 = -0.1

    sbend = converter.make_sbend(ns)
    assert sbend.id == ns.NAME
    assert sbend.l == ns.L
    assert sbend.angle == ns.ANGLE
    assert sbend.e1 == ns.E1
    assert sbend.e2 == ns.E2

@pytest.mark.skipif(MISSING_TFS_DEPENDENCY, reason=REASON)
def test_MADXLatticeConverter_make_rbend(converter, ns):
    ns.ANGLE = 0.5
    ns.TILT = 0.25
    ns.K1L = 1.3
    ns.K2L = 14.
    ns.FINT = 13.
    ns.FINTX = 12.
    ns.E1 = +0.1
    ns.E2 = -0.1


    rbend = converter.make_rbend(ns)
    assert rbend.id == ns.NAME
    assert rbend.l == ns.L
    assert rbend.angle == ns.ANGLE
    assert rbend.k1 == ns.K1L / ns.L
    assert rbend.k2 == ns.K2L / ns.L
    assert rbend.tilt == ns.TILT
    assert rbend.e1 == ns.E1 + ns.ANGLE / 2
    assert rbend.e2 == ns.E2 + ns.ANGLE / 2
    assert rbend.fint == ns.FINT
    assert rbend.fintx == ns.FINTX

@pytest.mark.skipif(MISSING_TFS_DEPENDENCY, reason=REASON)
def test_MADXLatticeConverter_make_quadrupole(converter, ns):
    ns.TILT = 0.25
    ns.K1L = 1.3

    quadrupole = converter.make_quadrupole(ns)
    assert quadrupole.id == ns.NAME
    assert quadrupole.l == ns.L
    assert quadrupole.k1 == ns.K1L / ns.L
    assert quadrupole.tilt == ns.TILT

@pytest.mark.skipif(MISSING_TFS_DEPENDENCY, reason=REASON)
def test_MADXLatticeConverter_make_sextupole(converter, ns):
    ns.TILT = 0.25
    ns.K2L = 1.3

    sextupole = converter.make_sextupole(ns)
    assert sextupole.id == ns.NAME
    assert sextupole.l == ns.L
    assert sextupole.k2 == ns.K2L / ns.L
    assert sextupole.tilt == ns.TILT

@pytest.mark.skipif(MISSING_TFS_DEPENDENCY, reason=REASON)
def test_MADXLatticeConverter_make_octupole(converter, ns):
    ns.TILT = 0.25
    ns.K3L = 1.3

    octupole = converter.make_octupole(ns)
    assert octupole.id == ns.NAME
    assert octupole.l == ns.L
    assert octupole.k3 == ns.K3L / ns.L
    assert octupole.tilt == ns.TILT

@pytest.mark.skipif(MISSING_TFS_DEPENDENCY, reason=REASON)
def test_MADXLatticeConverter_make_vkicker(converter, ns):
    ns.TILT = 0.25
    ns.VKICK = 1.3

    vcor = converter.make_vkicker(ns)
    assert vcor.id == ns.NAME
    assert vcor.l == ns.L
    assert vcor.angle == ns.VKICK

@pytest.mark.skipif(MISSING_TFS_DEPENDENCY, reason=REASON)
def test_MADXLatticeConverter_make_hkicker(converter, ns):
    ns.TILT = 0.25
    ns.HKICK = 1.3

    hcor = converter.make_hkicker(ns)
    assert hcor.id == ns.NAME
    assert hcor.l == ns.L
    assert hcor.angle == ns.HKICK
