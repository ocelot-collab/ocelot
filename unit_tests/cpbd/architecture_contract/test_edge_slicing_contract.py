import pytest

from ocelot.cpbd.elements import Cavity, Drift, SBend
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.transformations.cavity import CavityTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.transformations.transformation import TMTypes


def _tm_types(tm_list):
    return [tm.tm_type for tm in tm_list]


class BrokenFirstOrderEdgeAtom(Element):
    def __init__(self):
        super().__init__(eid="BROKEN_EDGE_FO", has_edge=True)
        self.l = 0.5


class BrokenFirstOrderEdgeElement(OpticElement):
    default_tm = TransferMap
    supported_tms = {TransferMap}

    def __init__(self):
        super().__init__(BrokenFirstOrderEdgeAtom(), tm=TransferMap, default_tm=TransferMap)


class FirstOrderOnlyEdgeAtom(Element):
    def __init__(self):
        super().__init__(eid="BROKEN_EDGE_SO", has_edge=True)
        self.l = 0.5

    def create_first_order_entrance_params(self, energy: float, delta_length: float = 0.0):
        return self.create_first_order_main_params(energy, delta_length)

    def create_first_order_exit_params(self, energy: float, delta_length: float = 0.0):
        return self.create_first_order_main_params(energy, delta_length)


class BrokenSecondOrderEdgeElement(OpticElement):
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM}

    def __init__(self):
        super().__init__(FirstOrderOnlyEdgeAtom(), tm=TransferMap, default_tm=TransferMap)


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


def test_sbend_ignore_edges_returns_only_main_slice():
    bend = SBend(l=1.2, angle=0.2, e1=0.05, e2=0.04, eid="B5")

    section = bend.get_section_tms(start_l=0.0, delta_l=bend.l, ignore_edges=True)

    assert _tm_types(section) == [TMTypes.MAIN]


def test_cavity_full_section_keeps_edge_sequence_with_cavity_tm():
    cavity = Cavity(l=0.7, freq=1.3e9, phi=20.0, v=0.005, eid="C_EDGE")

    section = cavity.get_section_tms(start_l=0.0, delta_l=cavity.l)

    assert _tm_types(section) == [TMTypes.ENTRANCE, TMTypes.MAIN, TMTypes.EXIT]
    assert all(isinstance(tm, CavityTM) for tm in section)


def test_has_edge_true_requires_first_order_edge_hooks_for_optics_path():
    with pytest.raises(RuntimeError, match="has has_edge=True"):
        BrokenFirstOrderEdgeElement()


def test_declared_edge_tm_without_edge_hooks_fails_clearly():
    element = BrokenSecondOrderEdgeElement()

    with pytest.raises(RuntimeError, match="has has_edge=True"):
        element.set_tm(SecondTM)
