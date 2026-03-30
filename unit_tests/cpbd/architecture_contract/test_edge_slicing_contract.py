import pytest

from ocelot.cpbd.elements import Drift, SBend
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.transformations.transformation import TMTypes


def _tm_types(tm_list):
    return [tm.tm_type for tm in tm_list]


def test_drift_section_uses_single_main_map():
    drift = Drift(l=1.5, eid="D1")

    section = drift.get_section_tms(start_l=0.2, delta_l=0.4)

    assert _tm_types(section) == [TMTypes.MAIN]
    assert section[0].delta_length == pytest.approx(0.4)


def test_sbend_full_section_includes_entrance_main_and_exit():
    bend = SBend(l=1.2, angle=0.2, e1=0.05, e2=0.04, eid="B1")

    section = bend.get_section_tms(start_l=0.0, delta_l=bend.l)

    assert _tm_types(section) == [TMTypes.ENTRANCE, TMTypes.MAIN, TMTypes.EXIT]
    assert section[0] is not bend.get_tm(TMTypes.ENTRANCE)
    assert section[1] is not bend.get_tm(TMTypes.MAIN)
    assert section[2] is not bend.get_tm(TMTypes.EXIT)


def test_sbend_leading_half_includes_entrance_and_main_only():
    bend = SBend(l=1.2, angle=0.2, e1=0.05, e2=0.04, eid="B2")

    section = bend.get_section_tms(start_l=0.0, delta_l=bend.l / 2.0)

    assert _tm_types(section) == [TMTypes.ENTRANCE, TMTypes.MAIN]
    assert section[1].delta_length == pytest.approx(bend.l / 2.0)


def test_sbend_trailing_half_includes_main_and_exit_only():
    bend = SBend(l=1.2, angle=0.2, e1=0.05, e2=0.04, eid="B3")

    section = bend.get_section_tms(start_l=bend.l / 2.0, delta_l=bend.l / 2.0)

    assert _tm_types(section) == [TMTypes.MAIN, TMTypes.EXIT]
    assert section[0].delta_length == pytest.approx(bend.l / 2.0)


def test_sbend_first_order_only_slice_uses_first_order_maps():
    bend = SBend(l=1.2, angle=0.2, e1=0.05, e2=0.04, tm=SecondTM, eid="B4")

    section = bend.get_section_tms(start_l=0.0, delta_l=bend.l / 2.0, first_order_only=True)

    assert _tm_types(section) == [TMTypes.ENTRANCE, TMTypes.MAIN]
    assert all(isinstance(tm, TransferMap) for tm in section)
