from functools import partial

from ocelot.cpbd.elements.xyquadrupole_atom import XYQuadrupoleAtom
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.transformations.transformation import TMTypes


class XYQuadrupole(OpticElement):
    """
    Linear transport about an off-axis quadrupole reference trajectory.

    ``x_offs`` and ``y_offs`` are the reference's entrance positions [m]
    relative to the magnetic axis, in the unrolled magnet frame. ``x_angle``
    and ``y_angle`` are its entrance tangent components dX/ds and dY/ds
    (the paper's theta_x0 and theta_y0, approximately angles in radians).
    They default to zero for parallel entry. ``k1`` [1/m^2] is positive for horizontal
    focusing, ``l`` [m] is reference arc length, and ``tilt`` [rad] rolls the
    complete map. Particle coordinates are deviations from the moving
    reference, so the zero particle stays zero even with nonzero offsets.
    The usual ``dx`` / ``dy`` alignment translations are separate parameters.

    Implements the analytical linear model in S. Tomin, "Beam Transport
    Through a Quadrupole with Offsets in Both Transverse Planes", revised v4,
    Eqs. (37)-(42) and Appendices A-B, with Ocelot's ``tau = -sigma`` convention.
    Curvatures vary along the reference orbit; curvature-squared focusing,
    fringe fields and higher-order particle transport are omitted. The model
    requires small reference angles and ``abs(k1) * X0(s)**2 << 1`` and
    ``abs(k1) * Y0(s)**2 << 1`` throughout the element. Map boundaries are
    planes normal to the reference, not fixed laboratory magnet faces.
    """
    default_tm = TransferMap
    supported_tms = {TransferMap}

    def __init__(self, l=0., x_offs=0.0, y_offs=0.0, k1=0.0, tilt=0.0, eid=None, tm=None,
                 x_angle=0.0, y_angle=0.0, **kwargs):
        super().__init__(XYQuadrupoleAtom(l=l, x_offs=x_offs, y_offs=y_offs, k1=k1, tilt=tilt, eid=eid,
                                       x_angle=x_angle, y_angle=y_angle, **kwargs), tm=tm)

    def get_section_tms(self, delta_l, start_l=0.0, ignore_edges=False, first_order_only=False):
        """Transport a slice along the original, continuously curved reference.

        The curvature depends on position in the magnet. A slice therefore
        uses ``R(start_l + delta_l) @ inv(R(start_l))``, not the entrance map
        of a shorter quadrupole with reset entrance conditions.
        """
        # Navigator obtains these distances by subtracting accumulated lattice
        # positions. Match its boundary tolerance before checking the interval.
        tolerance = 1e-10
        if start_l < -tolerance or start_l > self.l + tolerance or delta_l < -tolerance:
            raise ValueError("XYQuadrupole slice requires 0 <= start_l <= l and delta_l >= 0")
        start_l = min(max(start_l, 0.0), self.l)
        delta_l = min(max(delta_l, 0.0), self.l - start_l)
        return [TransferMap(
            create_tm_param_func=partial(
                self.element.create_first_order_main_params, start_length=start_l),
            delta_e_func=self.element.create_delta_e,
            tm_type=TMTypes.MAIN, length=self.l, delta_length=delta_l,
        )]
