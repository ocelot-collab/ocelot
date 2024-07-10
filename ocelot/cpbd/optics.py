__author__ = 'Sergey'


from numpy.linalg import inv

from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.beam import Twiss, twiss_iterable_to_df

from ocelot.cpbd.r_matrix import *
from ocelot.cpbd.tm_utils import SecondOrderMult
from ocelot.cpbd.transformations.second_order import SecondTM

_logger = logging.getLogger(__name__)


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
    """
    obj_list = []
    i = 0
    elem = lattice.sequence[i]
    L = elem.l
    obj_elem = obj0
    for z in z_array:
        while z > L:
            for tm in lattice.sequence[i].first_order_tms:
                obj_elem = tm * obj_elem
            i += 1
            elem = lattice.sequence[i]
            L += elem.l

        delta_l = z - (L - elem.l)
        first_order_tms = elem.get_section_tms(start_l=0.0, delta_l=delta_l, first_order_only=True)

        obj_z = obj_elem
        for tm in first_order_tms:
            obj_z = tm * obj_z

        obj_list.append(obj_z)
    return obj_list


def trace_obj(lattice, obj, nPoints=None, attach2elem=False):
    """
    track object through the lattice
    obj must be Twiss or Particle
    """

    if nPoints is None:
        obj_list = [obj]
        for e in lattice.sequence:
            for tm in e.first_order_tms:
                obj = tm * obj
                obj.id = e.id
                obj_list.append(obj)
            if attach2elem:
                e.tws = obj
    else:
        z_array = np.linspace(0, lattice.totalLen, nPoints, endpoint=True)
        obj_list = trace_z(lattice, obj, z_array)
    return obj_list


def periodic_twiss(tws, R):
    """
    initial conditions for a periodic Twiss solution
    """
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

    if abs(cosmx) >= 1 or abs(cosmy) >= 1:
        _logger.warning(" ************ periodic solution does not exist. return None ***********")
        return None
    sinmx = np.sign(R[0, 1]) * np.sqrt(1. - cosmx * cosmx)
    sinmy = np.sign(R[2, 3]) * np.sqrt(1. - cosmy * cosmy)

    tws.beta_x = abs(R[0, 1] / sinmx)
    tws.beta_y = abs(R[2, 3] / sinmy)

    tws.alpha_x = (R[0, 0] - R[1, 1]) / (2. * sinmx)  # X[0,0]

    tws.gamma_x = (1. + tws.alpha_x * tws.alpha_x) / tws.beta_x  # X[1,0]

    tws.alpha_y = (R[2, 2] - R[3, 3]) / (2 * sinmy)  # Y[0,0]
    tws.gamma_y = (1. + tws.alpha_y * tws.alpha_y) / tws.beta_y  # Y[1,0]

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


def twiss(lattice, tws0=None, nPoints=None, return_df=False, attach2elem=False):
    """
    twiss parameters calculation

    :param attach2elem: if True and nPoints=None Twiss will be attached to 'elem.tws', Twiss corresponds to the end of the element,
                        not recommended for standard use, but may be handy for small scripts
    :param return_df:
    :param lattice: lattice, MagneticLattice() object
    :param tws0: initial twiss parameters, Twiss() object. If None, function tries to find periodic solution.
    :param nPoints: number of points per cell. If None, then twiss parameters are calculated at the end of each element.
    :return: list of Twiss() objects
    """
    if tws0 is None:
        tws0 = lattice.periodic_twiss(tws0)

    if tws0.__class__ == Twiss:
        if tws0.beta_x == 0 or tws0.beta_y == 0:
            tws0 = lattice.periodic_twiss(tws0)
            if tws0 is None:
                _logger.info(' twiss: Twiss: no periodic solution')
                return None
        else:
            tws0.gamma_x = (1. + tws0.alpha_x ** 2) / tws0.beta_x
            tws0.gamma_y = (1. + tws0.alpha_y ** 2) / tws0.beta_y

        twiss_list = trace_obj(lattice, tws0, nPoints, attach2elem)

        if return_df:
            twiss_list = twiss_iterable_to_df(twiss_list)

        return twiss_list
    else:
        _logger.warning(' Twiss: no periodic solution. return None')
        return None


def twiss_fast(lattice, tws0=None):
    """
    twiss parameters calculation

    :param lattice: lattice, MagneticLattice() object
    :param tws0: initial twiss parameters, Twiss() object. If None, try to find periodic solution.
    :param nPoints: number of points per cell. If None, then twiss parameters are calculated at the end of each element.
    :return: list of Twiss() objects
    """
    if tws0 is None:
        tws0 = lattice.periodic_twiss(tws0)
    if tws0.__class__ == Twiss:
        if tws0.beta_x == 0 or tws0.beta_y == 0:
            tws0 = lattice.periodic_twiss(tws0)
            if tws0 is None:
                _logger.warning(' twiss_fast: Twiss: no periodic solution')
                return None
        else:
            tws0.gamma_x = (1. + tws0.alpha_x ** 2) / tws0.beta_x
            tws0.gamma_y = (1. + tws0.alpha_y ** 2) / tws0.beta_y

        obj_list = [tws0]
        for e in lattice.fast_seq:
            e.transfer_map.R = lambda x: e.transfer_map._r
            tws0 = e.transfer_map * tws0
            tws0.id = e.id
            obj_list.append(tws0)
        return obj_list
    else:
        _logger.warning(' twiss_fast: Twiss: no periodic solution')
        return None


def merge_maps(t_maps):
    tm0 = TransferMap()
    t_maps_new = []
    for tm in t_maps:
        if tm.__class__ == TransferMap:
            tm0 = tm * tm0
        else:
            t_maps_new.append(tm0)
            t_maps_new.append(tm)
            tm0 = TransferMap()
    t_maps_new.append(tm0)
    return t_maps_new


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
