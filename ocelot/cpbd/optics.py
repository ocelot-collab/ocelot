__author__ = 'Sergey'


from copy import deepcopy
from numpy.linalg import inv
import pandas as pd
from typing import Iterable
import numpy as np
import logging

from ocelot.cpbd.transformations.transfer_map import TransferMap

from ocelot.common.globals import m_e_GeV
from ocelot.cpbd.tm_utils import SecondOrderMult
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.beam import Twiss, twiss_iterable_to_df

logger = logging.getLogger(__name__)

class UnstableLatticeError(RuntimeError):
    """Raised when a lattice has no strictly stable periodic Twiss solution."""

    def __init__(self, cos_mu_x, cos_mu_y):
        self.cos_mu_x = cos_mu_x
        self.cos_mu_y = cos_mu_y
        unstable_planes = []
        for plane, cos_mu in (("x", cos_mu_x), ("y", cos_mu_y)):
            if not np.isfinite(cos_mu) or abs(cos_mu) >= 1:
                unstable_planes.append(f"{plane}: |Tr(R)/2|={abs(cos_mu):.6g}")
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
    """Return periodic initial Twiss parameters for a one-turn matrix."""
    tws = Twiss(tws)

    if R[5, 5] != 1:
        if tws.E == 0:
            raise TypeError("Lattice is contained Cavity. Argument 'tws' must be Twiss class with non zero energy 'tws.E'")

        g0 = tws.E / m_e_GeV
        g1 = np.sqrt(g0 ** 2 - 1 + R[5, 5] ** 2) / R[5, 5]
        k = np.sqrt(g1 / g0)
        R[0, 0] = R[0, 0] * k
        R[0, 1] = R[0, 1] * k
        R[1, 0] = R[1, 0] * k
        R[1, 1] = R[1, 1] * k
        R[2, 2] = R[2, 2] * k
        R[2, 3] = R[2, 3] * k
        R[3, 2] = R[3, 2] * k
        R[3, 3] = R[3, 3] * k

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

    tws.alpha_x = (R[0, 0] - R[1, 1]) / (2. * sinmx)  # X[0,0]


    tws.alpha_y = (R[2, 2] - R[3, 3]) / (2 * sinmy)  # Y[0,0]

    Hx = np.array([[R[0, 0] - 1, R[0, 1]], [R[1, 0], R[1, 1] - 1]])
    Hhx = np.array([[R[0, 5]], [R[1, 5]]])
    hh = np.dot(inv(-Hx), Hhx)
    tws.Dx = hh[0, 0]
    tws.Dxp = hh[1, 0]
    Hy = np.array([[R[2, 2] - 1, R[2, 3]], [R[3, 2], R[3, 3] - 1]])
    Hhy = np.array([[R[2, 5]], [R[3, 5]]])
    hhy = np.dot(inv(-Hy), Hhy)
    tws.Dy = hhy[0, 0]
    tws.Dyp = hhy[1, 0]
    return tws


def _validate_twiss_seed(tws0):
    if not isinstance(tws0, Twiss):
        raise TypeError(f"tws0 must be a Twiss instance, got {type(tws0).__name__}")
    if tws0.beta_x <= 0 or tws0.beta_y <= 0:
        raise ValueError(
            "tws0.beta_x and tws0.beta_y must be positive; "
            "use periodic_twiss(lattice, ...) to calculate periodic optics"
        )


def twiss(lattice, tws0, nPoints=None, return_df=False, attach2elem=False):
    """
    twiss parameters calculation

    :param attach2elem: if True and nPoints=None, attach the exit Twiss to every
                        element as ``elem.tws``. An iterable attaches only to
                        the selected elements. Every selected element instance
                        must identify exactly one lattice occurrence. This is
                        not recommended for standard use, but may be handy for
                        small scripts.
    :param return_df:
    :param lattice: lattice, MagneticLattice() object
    :param tws0: explicit initial Twiss parameters to propagate. ``beta_x``
                 and ``beta_y`` must be positive. Use :func:`periodic_twiss`
                 when the initial optics should be calculated from the lattice.
    :param nPoints: number of points per cell. If None, then twiss parameters are calculated at the end of each element.
    :return: list of Twiss() objects
    """
    attachment_targets = _resolve_attachment_targets(lattice, attach2elem)
    if attachment_targets and nPoints is not None:
        raise ValueError("attach2elem requires nPoints=None so element exits are available")

    _validate_twiss_seed(tws0)

    attachment_elements = [element for element, _index in attachment_targets]
    twiss_list = trace_obj(lattice, tws0, nPoints, attachment_elements)

    if return_df:
        twiss_list = twiss_iterable_to_df(twiss_list)

    return twiss_list


def periodic_twiss(lattice, tws0=None, nPoints=None, return_df=False, attach2elem=False):
    """Calculate and propagate periodic Twiss parameters through a lattice.

    ``tws0`` may provide beam energy, emittance, and other seed values. Its
    transverse Twiss parameters are replaced by the periodic solution.

    Raises
    ------
    UnstableLatticeError
        If either transverse plane has no strictly stable periodic solution.
    """

    periodic_seed = lattice.periodic_twiss(tws=tws0)
    return twiss(
        lattice,
        periodic_seed,
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
    _validate_twiss_seed(tws0)

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
        if not self.params.get('nkick') != self.nkick:
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
