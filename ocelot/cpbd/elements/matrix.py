from typing import List

import numpy as np

from ocelot.cpbd.elements.drift_atom import DriftAtom
from ocelot.cpbd.elements.matrix_atom import MatrixAtom
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transformation import Transformation
from ocelot.cpbd.transformations.transfer_map import TransferMap


class _MatrixApproxSliceAtom(DriftAtom):
    """
    Approximate optics-only slice used for partial Matrix sections.

    ``Matrix`` stores an exact full-element black-box map. For a partial
    section we intentionally do not try to factor the stored ``R`` / ``B`` /
    ``T`` into a physically meaningful interior slice. Instead we reuse the
    existing ``DriftAtom`` first-order path to create a simple optics
    interpolation map while preserving the wrapper-visible slice length,
    offsets/tilt, and linearly scaled ``delta_e``.
    """

    def __init__(self, source: MatrixAtom, l: float, delta_e: float):
        super().__init__(l=l, eid=source.id)
        self.delta_e = delta_e
        self.dx = source.dx
        self.dy = source.dy
        self.tilt = source.tilt

    def create_delta_e(self, total_length, delta_length=0.0):
        return self.delta_e


class Matrix(OpticElement):
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM}

    def __init__(self, l=0., delta_e=0, eid=None, tm=None, **kwargs):
        super().__init__(MatrixAtom(l=l, delta_e=delta_e, eid=eid, **kwargs), tm=tm)

    def get_section_tms(self, delta_l: float, start_l: float = 0.0, ignore_edges=False, first_order_only=False) -> List[Transformation]:
        total_length = self.element.l
        if start_l + delta_l > total_length or np.isclose(start_l + delta_l, total_length):
            delta_l = max(0.0, total_length - start_l)

        if start_l < 1e-10 and np.isclose(delta_l, total_length):
            return super().get_section_tms(delta_l=delta_l, start_l=start_l, ignore_edges=ignore_edges, first_order_only=first_order_only)

        if not first_order_only:
            # Active tracking slices are forbidden because Matrix does not
            # define a trustworthy interior nonlinear map. Only the full stored
            # black-box map is considered exact.
            raise RuntimeError(
                "Matrix supports partial slicing only for first_order_only=True optics interpolation; "
                "active tracking slices are not supported."
            )

        slice_delta_e = self.delta_e * delta_l / total_length if total_length != 0 else self.delta_e
        # For optics/Twiss sampling, approximate the interior of Matrix as a
        # drift-like first-order slice with scaled reference-energy gain.
        slice_atom = _MatrixApproxSliceAtom(source=self.element, l=delta_l, delta_e=slice_delta_e)
        return [TransferMap.from_element(slice_atom)]
