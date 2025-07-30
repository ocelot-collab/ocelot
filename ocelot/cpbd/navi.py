from ocelot.cpbd.physics_proc import RectAperture, EllipticalAperture
from ocelot.common.ocelog import *
from ocelot.cpbd.elements.aperture import Aperture
import numpy as np
from copy import deepcopy
_logger_navi = logging.getLogger(__name__ + ".navi")


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
    Navigator defines the step size (dz) for tracking and specifies which physical processes will be applied during each step.

    Attributes:
        lattice (MagneticLattice): The magnetic lattice to which the navigator is applied.
        unit_step (float): Unit step size for all physics processes, default is 1 meter.

    Methods:
        add_physics_proc(physics_proc, elem1, elem2):
            Adds a physical process to be applied between two elements.

            Parameters:
                physics_proc: The physical process to be applied (e.g., CSR, SpaceCharge, Wake, etc.).
                elem1: The element where the physical process begins (applied at the start of elem1).
                elem2: The element where the physical process ends (applied up to the start of elem2).

            Note:
                The physical process starts at the beginning of elem1 and ends at the beginning of elem2.
    """

    def __init__(self, lattice, unit_step=1):

        self.lat = lattice
        if len(self.lat.sequence) < 2:
            _logger_navi.warning("The lattice contains only one element. Adding a physics process might not work.")
        self.process_table = ProcessTable(self.lat)
        self.ref_process_table = deepcopy(self.process_table)

        self.z0 = 0.  # current position of navigator
        self.n_elem = 0  # current index of the element in lattice
        self.sum_lengths = 0.  # sum_lengths = Sum[lat.sequence[i].l, {i, 0, n_elem-1}]
        self.unit_step = unit_step  # unit step for physics processes
        self.proc_kick_elems = []
        self.kill_process = False  # for case when calculations are needed to terminated e.g. from gui
        self.inactive_processes = [] # processes are sometimes deactivated during tracking

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
        self.process_table.lat = self.lat
        self.inactive_processes = []
        self._update_references()
        for i, physics_proc in enumerate(self.process_table.proc_list):
            physics_proc.prepare(self.lat)


    def go_to_start(self):
        self.reset_position()

    def get_phys_procs(self):
        """
        method return list of all physics processes which were added

        :return: list, list of PhysProc(s)
        """
        return self.process_table.proc_list

    def add_physics_proc(self, physics_proc, elem1, elem2) -> None:
        """
        Adds a physical process to be applied between two elements in the lattice.

        :param physics_proc: The physical process to be applied, e.g., SpaceCharge, CSR, Wake, etc.
        :param elem1: The element in the lattice where the physical process starts (applied at the start of elem1).
        :param elem2: The element in the lattice where the physical process stops (applied up to the start of elem2).
                      elem2 can be the same as elem1.
        :return: None
        """
        _logger_navi.debug(" add_physics_proc: phys proc: " + physics_proc.__class__.__name__)
        self.process_table.add_physics_proc(physics_proc, elem1, elem2)
        self.ref_process_table = deepcopy(self.process_table)

    def add_physics_processes(self, processes, elem1s, elem2s) -> None:
        # Same as add physics_proc but much faster if you wanto add lots of processes.
        for process, e1, e2 in zip(processes, elem1s, elem2s):
            self.process_table.add_physics_proc(process, e1, e2)
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
            if elem.__class__ == Aperture:
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
        in case physics processes are applied and do not more needed they are
        removed from table.  They are moved to self.inactive_processes

        :param processes: list of processes are about to apply
        :return: None
        """
        for p in processes:
            if p in self.process_table.kick_proc_list:
                _logger_navi.debug(" Navigator.remove_used_processes: " + p.__class__.__name__)
                self.process_table.kick_proc_list.remove(p)
                self.process_table.proc_list.remove(p)
                self.inactive_processes.append(p)

    def jump_to(self, z: float) -> None:
        """
        Set the current position of the Navigator to some new z value.
        :param z: Position to jump to in the beamline.
        :return: None
        """
        self.go_to_start()
        if not z: # If jumping to 0.0 then we don't need to do anything.
            return
        if z < 0.0 or z >= self.lat.totalLen:
            raise ValueError(f"Invalid z value: {z}")


        self.z0 = z

        running_z = 0.0
        # update sum_lengths and n_elem for this new z point
        for i, element in enumerate(self.lat.sequence):
            running_z += element.l
            if running_z > self.z0:
                self.n_elem = i
                break
            self.sum_lengths += element.l

        # Remove any kick processes so we are consistent if we start
        # inside an element with a kick.
        for process in self.process_table.kick_proc_list:
            if process.s_stop < self.z0:
                self.process_table.kick_proc_list.remove(process)
                self.process_table.proc_list.remove(process)
                self.inactive_processes.append(process)

        # Update all process z0.
        for process in self.process_table.proc_list:
            process.z0 = z
        # update process table references to match those of the parent
        # Navigator class
        self._update_references()

    def get_next_step(self):
        while np.abs(self.z0 - self.lat.totalLen) > 1e-10:
            if self.kill_process:
                _logger_navi.info("Killing tracking ... ")
                return
            dz, proc_list, phys_steps = self.get_next()
            if self.z0 + dz > self.lat.totalLen:
                dz = self.lat.totalLen - self.z0

            yield self.get_map(dz), dz, proc_list, phys_steps

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

    def __str__(self):
        s = "Navigator: added physics processes: \n"
        for physproc in self.ref_process_table.proc_list:
            s += physproc.__class__.__name__ + " start: " + str(physproc.s_start) + "/ stop: " + str(physproc.s_stop) + "\n"
        return s

    def get_map(self, dz):
        nelems = len(self.lat.sequence)
        TM = []
        i = self.n_elem
        z1 = self.z0 + dz
        elem = self.lat.sequence[i]
        # navi.sum_lengths = np.sum([elem.l for elem in lattice.sequence[:i]])
        L = self.sum_lengths + elem.l
        while z1 + 1e-10 > L:

            dl = L - self.z0
            TM += elem.get_section_tms(start_l=self.z0 + elem.l - L, delta_l=dl)

            self.z0 = L
            dz -= dl
            if i >= nelems - 1:
                break

            i += 1
            elem = self.lat.sequence[i]
            L += elem.l

        if abs(dz) > 1e-10:
            TM += elem.get_section_tms(start_l=self.z0 + elem.l - L, delta_l=dz)

        self.z0 += dz
        self.sum_lengths = L - elem.l
        self.n_elem = i
        return TM

    def _update_references(self):
        # At initialisation the ProcessTable and Navi instances both
        # share the same lat instance.  However, ProcessTable is reset
        # (deepcopied) at various times and so we end up with the lats
        # having different addresses.  so we set some things here to
        # keep everything consistent.
        self.process_table.lat = self.lat
        for process in self.process_table.proc_list:
            indx0 = process.indx0
            indx1 = process.indx1
            process.start_elem = self.lat.sequence[indx0]
            process.end_elem = self.lat.sequence[indx1]
