__author__ = 'Sergey'

from copy import deepcopy
from ocelot.cpbd.transformations.multipole import MultipoleTM
from ocelot.cpbd.tm_params.second_order_params import SecondOrderParams

from numpy.linalg import inv

from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.beam import Twiss
from ocelot.cpbd.physics_proc import RectAperture, EllipticalAperture
from ocelot.cpbd.high_order import *

from ocelot.cpbd.r_matrix import *
from ocelot.cpbd.tm_utils import SecondOrderMult, transfer_maps_mult, unsym_matrix
from ocelot.cpbd.transformations.second_order import SecondTM

_logger = logging.getLogger(__name__)
_logger_navi = logging.getLogger(__name__ + ".navi")


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
        if not self.params.get('nKick') != self.nkick:
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

    Ra = np.eye(6)
    Ta = np.zeros((6, 6, 6))
    Ba = np.zeros((6, 1))
    E = energy
    for elem in lattice.sequence:
        for Rb, Bb, Tb, tm in zip(elem.R(E), elem.B(E), elem.T(E), elem.tms):
            Ra, Ta = transfer_maps_mult(Ra, Ta, Rb, Tb)
            Ba = np.dot(Rb, Ba) + Bb
            E += tm.get_delta_e()

    # TODO: Adding Attributes at runtime should be avoided
    lattice.E = E
    lattice.T_sym = Ta
    lattice.T = unsym_matrix(deepcopy(Ta))
    lattice.R = Ra
    lattice.B = Ba
    return Ra


def trace_z(lattice, obj0, z_array):
    """
    Z-dependent tracer (twiss(z) and particle(z))
    usage: twiss = trace_z(lattice,twiss_0, [1.23, 2.56, ...]) ,
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


def trace_obj(lattice, obj, nPoints=None):
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
    else:
        z_array = np.linspace(0, lattice.totalLen, nPoints, endpoint=True)
        obj_list = trace_z(lattice, obj, z_array)
    return obj_list


def periodic_twiss(tws, R):
    """
    initial conditions for a periodic Twiss solution
    """
    tws = Twiss(tws)

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


def twiss(lattice, tws0=None, nPoints=None):
    """
    twiss parameters calculation

    :param lattice: lattice, MagneticLattice() object
    :param tws0: initial twiss parameters, Twiss() object. If None, try to find periodic solution.
    :param nPoints: number of points per cell. If None, then twiss parameters are calculated at the end of each element.
    :return: list of Twiss() objects
    """
    if tws0 is None:
        tws0 = periodic_twiss(tws0, lattice_transfer_map(lattice, energy=0.))

    if tws0.__class__ == Twiss:
        if tws0.beta_x == 0 or tws0.beta_y == 0:
            R = lattice_transfer_map(lattice, tws0.E)
            tws0 = periodic_twiss(tws0, R)
            if tws0 is None:
                _logger.info(' twiss: Twiss: no periodic solution')
                return None
        else:
            tws0.gamma_x = (1. + tws0.alpha_x ** 2) / tws0.beta_x
            tws0.gamma_y = (1. + tws0.alpha_y ** 2) / tws0.beta_y

        twiss_list = trace_obj(lattice, tws0, nPoints)
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
        tws0 = periodic_twiss(tws0, lattice_transfer_map(lattice, energy=0.))
    if tws0.__class__ == Twiss:
        if tws0.beta_x == 0 or tws0.beta_y == 0:
            R = lattice_transfer_map(lattice, tws0.E)
            tws0 = periodic_twiss(tws0, R)
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


class ProcessTable:
    def __init__(self, lattice):
        self.proc_list = []
        self.kick_proc_list = []
        self.lat = lattice

    def searching_kick_proc(self, physics_proc, elem1):
        """
        function finds kick physics process. Kick physics process applies kick only once between two elements
        with zero length (e.g. Marker) or at the beginning of the element if it is the same element,
        others physics processes are applied during finite lengths.
        :return:
        """

        if (physics_proc.indx0 == physics_proc.indx1 or
                (physics_proc.indx0 + 1 == physics_proc.indx1 and elem1.l == 0)):
            physics_proc.indx1 = physics_proc.indx0
            physics_proc.s_stop = physics_proc.s_start
            self.kick_proc_list.append(physics_proc)
            if len(self.kick_proc_list) > 1:
                pos = np.array([proc.s_start for proc in self.kick_proc_list])
                indx = np.argsort(pos)
                self.kick_proc_list = [self.kick_proc_list[i] for i in indx]
        _logger_navi.debug(" searching_kick_proc: self.kick_proc_list.append(): " + str([p.__class__.__name__ for p in self.kick_proc_list]))

    def add_physics_proc(self, physics_proc, elem1, elem2):
        physics_proc.start_elem = elem1
        physics_proc.end_elem = elem2
        physics_proc.indx0 = self.lat.sequence.index(elem1)
        physics_proc.indx1 = self.lat.sequence.index(elem2)
        physics_proc.s_start = np.sum(np.array([elem.l for elem in self.lat.sequence[:physics_proc.indx0]]))
        physics_proc.s_stop = np.sum(np.array([elem.l for elem in self.lat.sequence[:physics_proc.indx1]]))
        self.searching_kick_proc(physics_proc, elem1)
        physics_proc.counter = physics_proc.step
        physics_proc.prepare(self.lat)

        _logger_navi.debug(" add_physics_proc: self.proc_list = " + str([p.__class__.__name__ for p in self.proc_list]) + ".append(" + physics_proc.__class__.__name__ + ")" +
                           "; start: " + str(physics_proc.indx0) + " stop: " + str(physics_proc.indx1))

        self.proc_list.append(physics_proc)


class Navigator:
    """
    Navigator defines step (dz) of tracking and which physical process will be applied during each step.
    lattice - MagneticLattice
    Attributes:
        unit_step = 1 [m] - unit step for all physics processes
    Methods:
        add_physics_proc(physics_proc, elem1, elem2)
            physics_proc - physics process, can be CSR, SpaceCharge or Wake,
            elem1 and elem2 - first and last elements between which the physics process will be applied.
    """

    def __init__(self, lattice):

        self.lat = lattice
        self.process_table = ProcessTable(self.lat)
        self.ref_process_table = deepcopy(self.process_table)

        self.z0 = 0.  # current position of navigator
        self.n_elem = 0  # current index of the element in lattice
        self.sum_lengths = 0.  # sum_lengths = Sum[lat.sequence[i].l, {i, 0, n_elem-1}]
        self.unit_step = 1  # unit step for physics processes
        self.proc_kick_elems = []
        self.kill_process = False  # for case when calculations are needed to terminated e.g. from gui

    def get_current_element(self):
        if self.n_elem < len(self.lat.sequence):
            return self.lat.sequence[self.n_elem]

    def reset_position(self):
        """
        method to reset Navigator position.
        :return:
        """
        _logger_navi.debug(" reset position")
        self.z0 = 0.  # current position of navigator
        self.n_elem = 0  # current index of the element in lattice
        self.sum_lengths = 0.  # sum_lengths = Sum[lat.sequence[i].l, {i, 0, n_elem-1}]
        self.process_table = deepcopy(self.ref_process_table)

    def go_to_start(self):
        self.reset_position()

    def get_phys_procs(self):
        """
        method return list of all physics processes which were added

        :return: list, list of PhysProc(s)
        """
        return self.process_table.proc_list

    def add_physics_proc(self, physics_proc, elem1, elem2):
        """
        Method adds Physics Process.

        :param physics_proc: PhysicsProc, e.g. SpaceCharge, CSR, Wake ...
        :param elem1: the element in the lattice where to start applying the physical process.
        :param elem2: the element in the lattice where to stop applying the physical process,
                        can be the same as starting element.
        :return:
        """
        _logger_navi.debug(" add_physics_proc: phys proc: " + physics_proc.__class__.__name__)
        self.process_table.add_physics_proc(physics_proc, elem1, elem2)
        self.ref_process_table = deepcopy(self.process_table)

    def activate_apertures(self, start=None, stop=None):
        """
        activate apertures if thea exist in the lattice from

        :param start: element,  activate apertures starting form element 'start' element
        :param stop: element, activate apertures up to 'stop' element
        :return:
        """
        id1 = self.lat.sequence.index(start) if start is not None else None
        id2 = self.lat.sequence.index(stop) if stop is not None else None
        for elem in self.lat.sequence[id1:id2]:
            # TODO: Move this logic to element class. __name__ is used temporary to break circular import.
            if elem.__class__.__name__ == "Aperture":
                if elem.type == "rect":
                    ap = RectAperture(xmin=-elem.xmax + elem.dx, xmax=elem.xmax + elem.dx,
                                      ymin=-elem.ymax + elem.dy, ymax=elem.ymax + elem.dy)
                    self.add_physics_proc(ap, elem, elem)
                elif elem.type == "ellipt":
                    ap = EllipticalAperture(xmax=elem.xmax, ymax=elem.ymax,
                                            dx=elem.dx, dy=elem.dy)
                    self.add_physics_proc(ap, elem, elem)

    def check_overjump(self, dz, processes, phys_steps):
        phys_steps_red = phys_steps - dz
        if len(processes) != 0:
            nearest_stop_elem = min([proc.indx1 for proc in processes])
            L_stop = np.sum(np.array([elem.l for elem in self.lat.sequence[:nearest_stop_elem]]))
            if self.z0 + dz > L_stop:
                dz = L_stop - self.z0

            # check if inside step dz there is another phys process

            proc_list_rest = [p for p in self.process_table.proc_list if p not in processes]
            start_pos_of_rest = np.array([proc.s_start for proc in proc_list_rest])
            start_pos = [pos for pos in start_pos_of_rest if self.z0 < pos < self.z0 + dz]
            if len(start_pos) > 0:
                start_pos.sort()
                dz = start_pos[0] - self.z0
                _logger_navi.debug(" check_overjump: there is phys proc inside step -> dz was decreased: dz = " + str(dz))

        phys_steps = phys_steps_red + dz

        # check kick processes
        kick_list = self.process_table.kick_proc_list
        kick_pos = np.array([proc.s_start for proc in kick_list])
        indx = np.argwhere(self.z0 < kick_pos)

        if 0 in kick_pos and self.z0 == 0 and self.n_elem == 0:
            indx0 = np.argwhere(self.z0 == kick_pos)
            indx = np.append(indx0, indx)
        if len(indx) != 0:
            kick_process = np.array([kick_list[i] for i in indx.flatten()])
            for i, proc in enumerate(kick_process):
                L_kick_stop = proc.s_start
                if self.z0 + dz > L_kick_stop:
                    dz = L_kick_stop - self.z0
                    phys_steps = phys_steps_red + dz
                    processes.append(proc)
                    phys_steps = np.append(phys_steps, 0)
                    continue
                elif self.z0 + dz == L_kick_stop:
                    processes.append(proc)
                    phys_steps = np.append(phys_steps, 0)
                else:
                    pass

        return dz, processes, phys_steps

    def get_proc_list(self):
        _logger_navi.debug(" get_proc_list: all phys proc = " + str([p.__class__.__name__ for p in self.process_table.proc_list]))
        proc_list = []
        for p in self.process_table.proc_list:
            if p.indx0 <= self.n_elem < p.indx1:
                proc_list.append(p)
        return proc_list

    def hard_edge_step(self, dz):
        # self.sum_lengths
        elem1 = self.lat.sequence[self.n_elem]
        L = self.sum_lengths + elem1.l
        if self.z0 + dz > L:
            dz = L - self.z0
        return dz

    def check_proc_bounds(self, dz, proc_list, phys_steps, active_process):
        for p in proc_list:
            if dz + self.z0 >= p.s_stop and p not in active_process:
                active_process.append(p)
                phys_steps = np.append(phys_steps, dz)

        return active_process, phys_steps

    def remove_used_processes(self, processes):
        """
        in case physics processes are applied and do not more needed they are removed from table

        :param processes: list of processes are about to apply
        :return: None
        """
        for p in processes:
            if p in self.process_table.kick_proc_list:
                _logger_navi.debug(" Navigator.remove_used_processes: " + p.__class__.__name__)
                self.process_table.kick_proc_list.remove(p)
                self.process_table.proc_list.remove(p)

    def get_next_step(self):
        while np.abs(self.z0 - self.lat.totalLen) > 1e-10:
            if self.kill_process:
                _logger.info("Killing tracking ... ")
                return
            dz, proc_list, phys_steps = self.get_next()
            if self.z0 + dz > self.lat.totalLen:
                dz = self.lat.totalLen - self.z0

            yield get_map(self.lat, dz, self), dz, proc_list, phys_steps

    def get_next(self):

        proc_list = self.get_proc_list()

        if len(proc_list) > 0:

            counters = np.array([p.counter for p in proc_list])
            step = counters.min()

            inxs = np.where(counters == step)

            processes = [proc_list[i] for i in inxs[0]]

            phys_steps = np.array([p.step for p in processes])*self.unit_step

            for p in proc_list:
                p.counter -= step
                if p.counter == 0:
                    p.counter = p.step

            dz = np.min(phys_steps)
        else:
            processes = proc_list
            n_elems = len(self.lat.sequence)
            if n_elems >= self.n_elem + 1:
                L = np.sum(np.array([elem.l for elem in self.lat.sequence[:self.n_elem + 1]]))
            else:
                L = self.lat.totalLen
            dz = L - self.z0
            phys_steps = np.array([])
        # check if dz overjumps the stop element
        dz, processes, phys_steps = self.check_overjump(dz, processes, phys_steps)
        processes, phys_steps = self.check_proc_bounds(dz, proc_list, phys_steps, processes)

        _logger_navi.debug(" Navigator.get_next: process: " + " ".join([proc.__class__.__name__ for proc in processes]))

        _logger_navi.debug(" Navigator.get_next: navi.z0=" + str(self.z0) + " navi.n_elem=" + str(self.n_elem) + " navi.sum_lengths="
                           + str(self.sum_lengths) + " dz=" + str(dz))

        _logger_navi.debug(" Navigator.get_next: element type=" + self.lat.sequence[self.n_elem].__class__.__name__ + " element name=" +
                           str(self.lat.sequence[self.n_elem].id))

        self.remove_used_processes(processes)

        return dz, processes, phys_steps


def get_map(lattice, dz, navi):
    nelems = len(lattice.sequence)
    TM = []
    i = navi.n_elem
    z1 = navi.z0 + dz
    elem = lattice.sequence[i]
    # navi.sum_lengths = np.sum([elem.l for elem in lattice.sequence[:i]])
    L = navi.sum_lengths + elem.l
    while z1 + 1e-10 > L:

        dl = L - navi.z0
        TM += elem.get_section_tms(start_l=navi.z0 + elem.l - L, delta_l=dl)

        navi.z0 = L
        dz -= dl
        if i >= nelems - 1:
            break

        i += 1
        elem = lattice.sequence[i]
        L += elem.l
        
    if abs(dz) > 1e-10:
        TM += elem.get_section_tms(start_l=navi.z0 + elem.l - L, delta_l=dz)

    navi.z0 += dz
    navi.sum_lengths = L - elem.l
    navi.n_elem = i
    return TM


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
