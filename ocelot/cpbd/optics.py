__author__ = 'Sergey'


from copy import deepcopy
import pandas as pd
from typing import Iterable, Literal
import numpy as np
import logging

from ocelot.cpbd.transformations.transfer_map import TransferMap

from ocelot.common.globals import m_e_GeV
from ocelot.cpbd.tm_utils import SecondOrderMult
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.beam import Particle, Twiss, twiss_iterable_to_df
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.magnetic_lattice import MagneticLattice

logger = logging.getLogger(__name__)

class UnstableLatticeError(RuntimeError):
    """Raised when a lattice has no strictly stable periodic Twiss solution."""

    def __init__(self, cos_mu_x: float, cos_mu_y: float, *, detail: str | None = None) -> None:
        self.cos_mu_x = cos_mu_x
        self.cos_mu_y = cos_mu_y
        unstable_planes = []
        for plane, cos_mu in (("x", cos_mu_x), ("y", cos_mu_y)):
            if not np.isfinite(cos_mu) or abs(cos_mu) >= 1:
                unstable_planes.append(f"{plane}: |Tr(R)/2|={abs(cos_mu):.6g}")
        if detail is not None:
            super().__init__(f"No strictly stable periodic Twiss solution exists ({detail}).")
            return
        details = ", ".join(unstable_planes)
        super().__init__(
            "No strictly stable periodic Twiss solution exists"
            f" ({details}; each value must be finite and less than 1)."
        )


def lattice_transfer_map(lattice, energy):
    """
    Function calculates transfer maps, the first and second orders (R, T), for the whole lattice.
    Second order matrices are attached to lattice object:
    lattice.T_sym - symmetric second order matrix
    lattice.T - second order matrix
    lattice.R - linear R matrix

    :param lattice: MagneticLattice
    :param energy: the initial electron beam energy [GeV]
    :return: R - matrix
    """
    Ba, Ra, Ta = lattice.transfer_maps(energy)

    # TODO: Adding Attributes at runtime should be avoided
    lattice.T_sym = Ta
    lattice.T = Ta #unsym_matrix(deepcopy(Ta))
    lattice.R = Ra
    lattice.B = Ba
    return Ra


def trace_z(lattice, obj0, z_array):
    """
    Z-dependent tracer (twiss(z) and particle(z))
    usage: twiss = trace_z(lattice, twiss_0, [1.23, 2.56, ...]) ,
    to calculate Twiss params at 1.23m, 2.56m etc.

    Positions are measured from the lattice entrance. Each sample is
    propagated independently from its element's entrance state, without
    modifying obj0.
    """
    eps = 1e-12
    obj_list = []
    i = 0
    elem = lattice.sequence[i]
    L = elem.l
    obj_elem = deepcopy(obj0)
    for z in z_array:
        while z > L + eps and i + 1 < len(lattice.sequence):
            for tm in lattice.sequence[i].first_order_tms:
                obj_elem = tm * obj_elem
            i += 1
            elem = lattice.sequence[i]
            L += elem.l

        delta_l = z - (L - elem.l)
        if delta_l < 0:   # safeguard against floating-point issues
            delta_l = 0.0
        elif delta_l > elem.l:
            delta_l = elem.l
        first_order_tms = elem.get_section_tms(start_l=0.0, delta_l=delta_l, first_order_only=True)

        # Particle transfer-map multiplication modifies its input before
        # returning a copy. Preserve the entrance state for later samples.
        obj_z = deepcopy(obj_elem)
        for tm in first_order_tms:
            obj_z = tm * obj_z

        obj_list.append(obj_z)
    return obj_list


def _resolve_attachment_targets(lattice, attach2elem):
    """Return unique ``(element, sequence_index)`` attachment targets."""

    if attach2elem is None or attach2elem is False:
        return []

    if attach2elem is True:
        requested = []
        seen = set()
        for element in lattice.sequence:
            identity = id(element)
            if identity not in seen:
                requested.append(element)
                seen.add(identity)
    else:
        if isinstance(attach2elem, (bool, np.bool_)):
            if not attach2elem:
                return []
            requested = list(lattice.sequence)
        else:
            try:
                requested = list(attach2elem)
            except TypeError as exc:
                raise TypeError(
                    "attach2elem must be a bool or an iterable of element instances"
                ) from exc

    targets = []
    seen = set()
    for element in requested:
        identity = id(element)
        if identity in seen:
            continue
        index = lattice.resolve_element_index(element)
        targets.append((element, index))
        seen.add(identity)
    return targets


def trace_obj(lattice, obj, nPoints=None, attach2elem=False):
    """
    track object through the lattice
    obj must be Twiss or Particle
    """

    attachment_targets = _resolve_attachment_targets(lattice, attach2elem)
    if attachment_targets and nPoints is not None:
        raise ValueError("attach2elem requires nPoints=None so element exits are available")

    if nPoints is None:
        obj_list = [obj]
        element_exit_objects = []
        for e in lattice.sequence:
            for tm in e.first_order_tms:
                obj = tm * obj
                obj.id = e.id
                obj_list.append(obj)
            element_exit_objects.append(obj)

        # Attach only after successful propagation so a failed calculation
        # cannot leave a partially updated lattice.
        for element, index in attachment_targets:
            element.tws = element_exit_objects[index]
    else:
        z_array = np.linspace(0, lattice.totalLen, nPoints, endpoint=True)
        obj_list = trace_z(lattice, obj, z_array)
    return obj_list


def _periodic_twiss_from_matrix(tws, R):
    """Match two uncoupled Courant–Snyder ellipses to the design map."""
    tws = Twiss(tws)
    roundoff = 64 * np.finfo(float).eps * max(1., np.linalg.norm(R[:4, :4]))
    if np.any(abs(R[:2, 2:4]) > roundoff) or np.any(abs(R[2:4, :2]) > roundoff):
        raise ValueError("Periodic design optics requires an uncoupled transverse one-turn map; "
                         "use periodic_twiss(..., reference_orbit='closed') for coupled rings")

    if R[5, 5] != 1:
        if tws.E == 0:
            raise TypeError("Lattice is contained Cavity. Argument 'tws' must be Twiss class with non zero energy 'tws.E'")

        g0 = tws.E / m_e_GeV
        g1 = np.sqrt(g0 ** 2 - 1 + R[5, 5] ** 2) / R[5, 5]
        k = np.sqrt(g1 / g0)
        R = R.copy()
        R[:4, :4] *= k

    cosmx = (R[0, 0] + R[1, 1]) / 2.
    cosmy = (R[2, 2] + R[3, 3]) / 2.

    if (
        not np.isfinite(cosmx)
        or not np.isfinite(cosmy)
        or abs(cosmx) >= 1
        or abs(cosmy) >= 1
    ):
        raise UnstableLatticeError(cosmx, cosmy)
    sinmx = np.sign(R[0, 1]) * np.sqrt(1. - cosmx * cosmx)
    sinmy = np.sign(R[2, 3]) * np.sqrt(1. - cosmy * cosmy)

    tws.beta_x = abs(R[0, 1] / sinmx)
    tws.beta_y = abs(R[2, 3] / sinmy)

    tws.alpha_x = (R[0, 0] - R[1, 1]) / (2. * sinmx)
    tws.alpha_y = (R[2, 2] - R[3, 3]) / (2. * sinmy)
    tws.xy = tws.xpy = tws.ypx = tws.pxpy = 0.
    tws._unit_moments = None
    tws._update_plane_moments()
    tws.Dx, tws.Dxp, tws.Dy, tws.Dyp = np.linalg.solve(np.eye(4) - R[:4, :4], R[:4, 5])
    return tws


def _prepare_entrance_twiss(tws0: Twiss, reference_orbit: Particle | None = None) -> Twiss:
    """Form the entrance optics sample while preserving the caller's beam moments."""
    if not isinstance(tws0, Twiss):
        raise TypeError(f"tws0 must be a Twiss instance, got {type(tws0).__name__}")
    if tws0.beta_x <= 0 or tws0.beta_y <= 0:
        raise ValueError(
            "tws0.beta_x and tws0.beta_y must be positive; "
            "use periodic_twiss(lattice, ...) to calculate periodic optics"
        )
    if reference_orbit is None:
        entrance_twiss = deepcopy(tws0)
    else:
        from ocelot.cpbd.orbit_optics import twiss_at_launch

        entrance_twiss = twiss_at_launch(tws0, reference_orbit)
    entrance_twiss._update_plane_moments()
    return entrance_twiss


def twiss(
    lattice: MagneticLattice,
    tws0: Twiss,
    nPoints: int | None = None,
    return_df: bool = False,
    attach2elem: bool | Iterable[OpticElement] | None = False,
    *,
    reference_orbit: Particle | None = None,
) -> list[Twiss] | pd.DataFrame:
    """Propagate entrance Twiss parameters through a beamline.

    By default, each element supplies its first-order design map. To calculate
    focusing along a steered beam, pass its entrance particle. Ocelot tracks
    that particle and transports the beam ellipse with the derivative of the
    selected tracking map at each point on its trajectory. For example, an
    off-axis sextupole contributes linear focusing when its tracking map
    includes the sextupole field.

    In an uncoupled line, ``mux`` and ``muy`` are the Courant–Snyder phase
    advances [rad]. To obtain normal-mode phases in a coupled ring, use
    :func:`periodic_twiss` with ``reference_orbit="closed"``.

    :param lattice: beamline in tracking order.
    :param tws0: entrance optics with positive ``beta_x`` and ``beta_y`` [m].
                 Supply both geometric emittances [m rad] for physical beam
                 sizes. In coupled transport their ratio affects projected beta.
    :param nPoints: number of equally spaced samples along the whole beamline.
                    The default includes the entrance and the end of each map.
    :param return_df: return a pandas DataFrame with one row per sample.
    :param attach2elem: with ``nPoints=None``, pass ``True`` or an iterable of
                        elements to store exit optics as ``elem.tws``. Each
                        selected object must occur exactly once in the lattice.
    :param reference_orbit: ``None`` for design optics. Pass ``Particle()`` for
                            a zero launch at ``tws0.E`` [GeV], or, for example,
                            ``Particle(x=1e-3)`` for a 1 mm horizontal launch.
    :return: list of Twiss samples, or a DataFrame. At every sample,
             ``xx = beta_x * emit_x`` and ``yy = beta_y * emit_y`` are the
             betatron beam variances [m²]. With a particle launch,
             ``x, xp, y, yp, tau, p`` give the tracked orbit.
    """
    attachment_targets = _resolve_attachment_targets(lattice, attach2elem)
    if attachment_targets and nPoints is not None:
        raise ValueError("attach2elem requires nPoints=None so element exits are available")

    entrance_twiss = _prepare_entrance_twiss(tws0, reference_orbit)

    if reference_orbit is None:
        attachment_elements = [element for element, _index in attachment_targets]
        twiss_list = trace_obj(lattice, entrance_twiss, nPoints, attachment_elements)
    else:
        from ocelot.cpbd.orbit_optics import trace_orbit_optics

        twiss_list = trace_orbit_optics(lattice, entrance_twiss, nPoints, attachment_targets)

    if return_df:
        twiss_list = twiss_iterable_to_df(twiss_list)

    return twiss_list


def periodic_twiss(
    lattice: MagneticLattice,
    tws0: Twiss | None = None,
    nPoints: int | None = None,
    return_df: bool = False,
    attach2elem: bool | Iterable[OpticElement] | None = False,
    *,
    reference_orbit: Particle | Literal["closed"] | None = None,
) -> list[Twiss] | pd.DataFrame:
    """Match the entrance optics to one turn and propagate it around the ring.

    The default uses the design maps and requires an uncoupled transverse
    one-turn matrix. For a coupled ring or a distorted orbit, pass
    ``reference_orbit="closed"``. The active tracking maps then determine the
    transverse closed orbit and the two normal modes of small oscillations
    about it. This calculation assumes fixed momentum and time-independent
    transverse motion. The arrival time may slip from turn to turn.

    :param lattice: one turn of the ring or one repeating cell.
    :param tws0: optional beam parameters giving the energy [GeV] and
                 geometric emittances [m rad]. In closed-orbit optics,
                 ``emit_x`` and ``emit_y`` weight the horizontal-like and vertical-like modes. The
                 returned beta, alpha and emittances describe the beam
                 projected onto the local horizontal and vertical planes.
    :param reference_orbit: ``None`` for design optics. Use ``"closed"`` to
                            start the orbit search at zero or pass a
                            ``Particle`` to set the transverse starting guess,
                            entrance ``tau`` [m] and ``p = (E-E0)/(p0*c)``.
    :return: Twiss samples along the ring. Sampling, DataFrame output and
             element attachment follow :func:`twiss`. In closed-orbit optics,
             ``mux`` and ``muy`` [rad] follow the phases of the horizontal and
             vertical position projections of the respective normal modes.

    Raises
    ------
    UnstableLatticeError
        If the transverse motion has no stable periodic solution.
    ValueError
        If the design one-turn matrix couples the planes, or the requested
        closed orbit or its symplectic transverse normal modes cannot be found.
    """

    if reference_orbit is not None:
        from ocelot.cpbd.orbit_optics import closed_orbit_optics, trace_orbit_optics

        targets = _resolve_attachment_targets(lattice, attach2elem)
        if targets and nPoints is not None:
            raise ValueError("attach2elem requires nPoints=None so element exits are available")
        matched_twiss, modes = closed_orbit_optics(lattice, tws0, reference_orbit)
        samples = trace_orbit_optics(lattice, matched_twiss, nPoints, targets, modes)
        return twiss_iterable_to_df(samples) if return_df else samples

    matched_twiss = lattice.periodic_twiss(tws=tws0)
    return twiss(
        lattice,
        matched_twiss,
        nPoints=nPoints,
        return_df=return_df,
        attach2elem=attach2elem,
    )


def twiss_fast(lattice, tws0):
    """
    twiss parameters calculation

    :param lattice: lattice, MagneticLattice() object
    :param tws0: explicit initial Twiss parameters.
    :param nPoints: number of points per cell. If None, then twiss parameters are calculated at the end of each element.
    :return: list of Twiss() objects
    """
    tws0 = _prepare_entrance_twiss(tws0)

    obj_list = [tws0]
    for e in lattice.fast_seq:
        e.transfer_map.R = lambda x: e.transfer_map._r
        tws0 = e.transfer_map * tws0
        tws0.id = e.id
        obj_list.append(tws0)
    return obj_list



class MethodTM:
    """
    The class creates a transfer map for elements that depend on user-defined parameters ("parameters").
    By default, the parameters = {"global": TransferMap}, which means that all elements will have linear transfer maps.
    You can also specify different transfer maps for any type of element.

    Example:
    --------
    # use linear matrices for all elements except Sextupole which will have nonlinear kick map (KickTM)
    method = MethodTM()
    method.global_method = TransferMap
    method.params[Sextupole] = KickTM

    # All elements are assigned matrices of the second order.
    # For elements for which there are no matrices of the second order are assigned default matrices, e.g. linear matrices.
    method2 = MethodTM()
    method2.global_method = SecondTM

    """

    def __init__(self, params=None):
        logger.warning("obsolete, use dictionary instead: {'global': SecondTM}")
        if params is None:
            self.params = {'global': TransferMap}
        else:
            self.params = params

        if "global" in self.params:
            self.global_method = self.params['global']
        else:
            self.global_method = TransferMap
        self.sec_order_mult = SecondOrderMult()
        self.nkick = self.params['nkick'] if 'nkick' in self.params else 1

    def to_dict(self):
        res = self.params
        if self.params.get('global') != self.global_method:
            res['global'] = self.global_method
        if self.params.get('nkick', 1) != self.nkick:
            res['nkick'] = self.nkick

        # OLD BEHAVIOR: old CorrectorTM has been splitted in First Order and Second Order to keep
        # the old behavior VCor's and Hcor's tm is set to SecondTM which is equal to
        # the old CorrectorTM.
        if not res.get('Vcor'):
            res['Vcor'] = SecondTM

        if not res.get('Hcor'):
            res['Hcor'] = SecondTM
        return res

'''
returns two solutions for a periodic fodo, given the mean beta
initial betas are at the center of the focusing quad
'''


def fodo_parameters(betaXmean=36.0, L=10.0, verbose=False):
    lquad = 0.001

    kap1 = np.sqrt(1.0 / 2.0 * (
        (betaXmean / L) * (betaXmean / L) + (betaXmean / L) * np.sqrt(-4.0 + (betaXmean / L) * (betaXmean / L))))
    kap2 = np.sqrt(1.0 / 2.0 * (
        (betaXmean / L) * (betaXmean / L) - (betaXmean / L) * np.sqrt(-4.0 + (betaXmean / L) * (betaXmean / L))))

    k = 1.0 / (lquad * L * kap2)

    f = 1.0 / (k * lquad)

    kappa = f / L
    betaMax = np.array(
        (L * kap1 * (kap1 + 1) / np.sqrt(kap1 * kap1 - 1), L * kap2 * (kap2 + 1) / np.sqrt(kap2 * kap2 - 1)))
    betaMin = np.array(
        (L * kap1 * (kap1 - 1) / np.sqrt(kap1 * kap1 - 1), L * kap2 * (kap2 - 1) / np.sqrt(kap2 * kap2 - 1)))
    betaMean = np.array(
        (L * kap2 * kap2 / (np.sqrt(kap2 * kap2 - 1.0)), L * kap1 * kap1 / (np.sqrt(kap1 * kap1 - 1.0))))
    k = np.array((1.0 / (lquad * L * kap1), 1.0 / (lquad * L * kap2)))

    if verbose:
        print('********* calculating fodo parameters *********')
        print('fodo parameters:')
        print('k*l=', k * lquad)
        print('f=', L * kap1, L * kap2)
        print('kap1=', kap1)
        print('kap2=', kap2)
        print('betaMax=', betaMax)
        print('betaMin=', betaMin)
        print('betaMean=', betaMean)
        print('*********                             *********')

    return k * lquad, betaMin, betaMax, betaMean
