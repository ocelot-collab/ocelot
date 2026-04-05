import numpy as np

from ocelot.cpbd.beam import Twiss
from ocelot.cpbd.elements import Cavity, Drift, Quadrupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.optics import twiss
from ocelot.cpbd.twiss_linear import (
    make_twiss_state,
    periodic_twiss_from_r,
    trace_linear_twiss,
    twiss_step_from_r,
)


def _state_from_twiss(tws):
    return make_twiss_state(
        beta_x=tws.beta_x,
        alpha_x=tws.alpha_x,
        beta_y=tws.beta_y,
        alpha_y=tws.alpha_y,
        Dx=tws.Dx,
        Dxp=tws.Dxp,
        Dy=tws.Dy,
        Dyp=tws.Dyp,
        mux=tws.mux,
        muy=tws.muy,
        E=tws.E,
        s=tws.s,
        p=tws.p,
    )


def _assert_state_matches_twiss(state, tws):
    np.testing.assert_allclose(state["beta_x"], tws.beta_x)
    np.testing.assert_allclose(state["alpha_x"], tws.alpha_x)
    np.testing.assert_allclose(state["beta_y"], tws.beta_y)
    np.testing.assert_allclose(state["alpha_y"], tws.alpha_y)
    np.testing.assert_allclose(state["Dx"], tws.Dx)
    np.testing.assert_allclose(state["Dxp"], tws.Dxp)
    np.testing.assert_allclose(state["Dy"], tws.Dy)
    np.testing.assert_allclose(state["Dyp"], tws.Dyp)
    np.testing.assert_allclose(state["mux"], tws.mux)
    np.testing.assert_allclose(state["muy"], tws.muy)
    np.testing.assert_allclose(state["E"], tws.E)
    np.testing.assert_allclose(state["s"], tws.s)


def test_twiss_step_from_r_matches_twiss_track():
    tws0 = Twiss(
        beta_x=11.0,
        alpha_x=1.3,
        beta_y=14.0,
        alpha_y=-0.4,
        Dx=0.2,
        Dxp=0.03,
        Dy=0.0,
        Dyp=0.01,
        E=0.0,
    )
    quad = Quadrupole(l=0.35, k1=1.4, tilt=0.12)
    R = quad.linear_r_full(energy=0.0)

    state = twiss_step_from_r(R, _state_from_twiss(tws0))
    ref = Twiss.track(R, tws0)

    _assert_state_matches_twiss(state, ref)


def test_trace_linear_twiss_matches_standard_twiss_at_element_boundaries():
    q1 = Quadrupole(l=0.25, k1=0.9, tilt=0.08, eid="Q1")
    d1 = Drift(l=0.6, eid="D1")
    q2 = Quadrupole(l=0.25, k1=-0.7, eid="Q2")
    lat = MagneticLattice([q1, d1, q2, d1])

    tws0 = Twiss(
        beta_x=12.0,
        alpha_x=0.5,
        beta_y=15.0,
        alpha_y=-0.2,
        Dx=0.1,
        Dxp=0.0,
        Dy=0.0,
        Dyp=0.0,
        E=0.0,
    )

    ref = twiss(lat, Twiss(tws0))
    states = trace_linear_twiss(lat, _state_from_twiss(tws0))

    assert len(states) == len(ref)
    for state, tws in zip(states, ref):
        _assert_state_matches_twiss(state, tws)


def test_trace_linear_twiss_supports_overrides_without_mutation():
    q1 = Quadrupole(l=0.25, k1=0.9, eid="Q1")
    d1 = Drift(l=0.6, eid="D1")
    lat = MagneticLattice([q1, d1])

    tws0 = Twiss(beta_x=12.0, alpha_x=0.5, beta_y=15.0, alpha_y=-0.2, E=0.0)
    state0 = _state_from_twiss(tws0)

    states = trace_linear_twiss(lat, state0, overrides={q1: {"k1": 1.25, "l": 0.35}})
    ref_lat = MagneticLattice([Quadrupole(l=0.35, k1=1.25, eid="Q1"), Drift(l=0.6, eid="D1")])
    ref = twiss(ref_lat, Twiss(tws0))

    assert q1.k1 == 0.9
    assert q1.l == 0.25
    _assert_state_matches_twiss(states[-1], ref[-1])


def test_periodic_twiss_from_r_matches_existing_periodic_twiss():
    d2 = Drift(l=0.5)
    qf = Quadrupole(l=0.2, k1=0.3)
    qdh = Quadrupole(l=0.1, k1=-0.3)
    lat = MagneticLattice((qdh, d2, Drift(l=0.5), d2, qf, d2, Drift(l=0.5), d2, qdh))

    energy = 0.005
    ref = lat.periodic_twiss(tws=Twiss(E=energy))
    R = lat.transfer_maps(energy=energy)[1]
    state = periodic_twiss_from_r(R, energy=energy)

    assert state["stable_x"]
    assert state["stable_y"]
    np.testing.assert_allclose(state["beta_x"], ref.beta_x)
    np.testing.assert_allclose(state["alpha_x"], ref.alpha_x)
    np.testing.assert_allclose(state["beta_y"], ref.beta_y)
    np.testing.assert_allclose(state["alpha_y"], ref.alpha_y)
    np.testing.assert_allclose(state["Dx"], ref.Dx)
    np.testing.assert_allclose(state["Dxp"], ref.Dxp)
    np.testing.assert_allclose(state["Dy"], ref.Dy)
    np.testing.assert_allclose(state["Dyp"], ref.Dyp)


def test_trace_linear_twiss_matches_standard_twiss_with_cavity_energy_gain():
    d1 = Drift(l=0.5, eid="D1")
    c1 = Cavity(l=0.35, v=0.02, phi=20.0, freq=1.3e9, eid="C1")
    q1 = Quadrupole(l=0.2, k1=0.7, eid="Q1")
    lat = MagneticLattice([d1, c1, q1, d1])

    tws0 = Twiss(
        beta_x=11.0,
        alpha_x=0.4,
        beta_y=13.0,
        alpha_y=-0.3,
        Dx=0.02,
        Dxp=0.0,
        Dy=0.0,
        Dyp=0.0,
        E=0.13,
    )

    ref = twiss(lat, Twiss(tws0))
    states = trace_linear_twiss(lat, _state_from_twiss(tws0))

    ref_boundary_states = [ref[0]]
    ref_index = 0
    for elem in lat.sequence:
        ref_index += len(elem.first_order_tms)
        ref_boundary_states.append(ref[ref_index])

    assert len(states) == len(ref_boundary_states)
    for state, tws in zip(states, ref_boundary_states):
        _assert_state_matches_twiss(state, tws)


def test_trace_linear_twiss_supports_cavity_overrides_without_mutation():
    c1 = Cavity(l=0.35, v=0.02, phi=20.0, freq=1.3e9, eid="C1")
    d1 = Drift(l=0.5, eid="D1")
    lat = MagneticLattice([c1, d1])

    tws0 = Twiss(beta_x=11.0, alpha_x=0.4, beta_y=13.0, alpha_y=-0.3, E=0.13)
    states = trace_linear_twiss(lat, _state_from_twiss(tws0), overrides={c1: {"v": 0.03, "phi": 15.0}})

    ref_lat = MagneticLattice([Cavity(l=0.35, v=0.03, phi=15.0, freq=1.3e9, eid="C1"), Drift(l=0.5, eid="D1")])
    ref = twiss(ref_lat, Twiss(tws0))

    assert c1.v == 0.02
    assert c1.phi == 20.0
    _assert_state_matches_twiss(states[-1], ref[-1])
