import numpy as np
import pytest

from ocelot.cpbd.elements import Cavity
from ocelot.cpbd.transformations.transformation import TMTypes


def _main_tm(tm_list):
    return next(tm for tm in tm_list if tm.tm_type == TMTypes.MAIN)


def test_cavity_only_main_map_changes_reference_energy():
    cavity = Cavity(l=1.2, v=0.03, phi=30.0, freq=1.3e9, eid="C3")
    expected_delta_e = cavity.v * np.cos(np.deg2rad(cavity.phi))

    delta_e_by_tm_type = {tm.tm_type: tm.get_delta_e() for tm in cavity.tms}

    assert delta_e_by_tm_type[TMTypes.ENTRANCE] == 0.0
    assert delta_e_by_tm_type[TMTypes.EXIT] == 0.0
    assert delta_e_by_tm_type[TMTypes.MAIN] == pytest.approx(expected_delta_e)


def test_cavity_slice_energy_gain_scales_with_slice_length():
    cavity = Cavity(l=1.2, v=0.03, phi=30.0, freq=1.3e9, eid="C4")
    expected_full_delta_e = cavity.v * np.cos(np.deg2rad(cavity.phi))

    leading_slice = cavity.get_section_tms(start_l=0.0, delta_l=cavity.l / 4.0)
    trailing_slice = cavity.get_section_tms(start_l=cavity.l / 4.0, delta_l=3.0 * cavity.l / 4.0)

    leading_main = _main_tm(leading_slice)
    trailing_main = _main_tm(trailing_slice)

    assert leading_main.get_delta_e() == pytest.approx(expected_full_delta_e / 4.0)
    assert trailing_main.get_delta_e() == pytest.approx(3.0 * expected_full_delta_e / 4.0)
    assert leading_main.get_delta_e() + trailing_main.get_delta_e() == pytest.approx(expected_full_delta_e)


def test_cavity_remove_coupler_kick_rebuilds_edge_maps():
    cavity = Cavity(
        l=0.346,
        v=0.0025,
        phi=180.0,
        freq=3.9e9,
        vx_up=0.1,
        vy_up=0.02,
        vxx_up=0.03,
        vxy_up=0.04,
        vx_down=0.05,
        vy_down=0.06,
        vxx_down=0.07,
        vxy_down=0.08,
        eid="C5",
    )
    energy = 0.13

    entrance_before = cavity.first_order_tms[0].get_params(energy)
    exit_before = cavity.first_order_tms[2].get_params(energy)
    entrance_b_before = entrance_before.B.copy()
    exit_b_before = exit_before.B.copy()
    entrance_r_before = entrance_before.R.copy()
    exit_r_before = exit_before.R.copy()

    cavity.remove_coupler_kick()

    entrance_after = cavity.first_order_tms[0].get_params(energy)
    exit_after = cavity.first_order_tms[2].get_params(energy)

    assert not np.allclose(entrance_b_before, entrance_after.B)
    assert not np.allclose(exit_b_before, exit_after.B)
    assert not np.allclose(entrance_r_before, entrance_after.R)
    assert not np.allclose(exit_r_before, exit_after.R)

    np.testing.assert_allclose(entrance_after.B, np.zeros((6, 1)))
    np.testing.assert_allclose(exit_after.B, np.zeros((6, 1)))
    np.testing.assert_allclose(entrance_after.R, np.eye(6))
    np.testing.assert_allclose(exit_after.R, np.eye(6))
