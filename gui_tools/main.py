"""
Tool for S2E simulation.
S.Tomin. XFEL/DESY. 2017
"""
import sys, os
path = os.path.realpath(__file__)
indx = path.find("gui_tools")
print("PATH", os.path.realpath(__file__), path[:indx])
sys.path.append(path[:indx])
from ocelot import *
from accelerator.s2e_sections import sections
from gui_tools.gui.gui_main import OcelotInterfaceWindow
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from PyQt5 import QtCore
from PyQt5.QtWidgets import QMainWindow, QApplication
import time
import numpy as np
from copy import deepcopy
from threading import Thread, Event
from ocelot.optimizer.mint.xfel_interface import *
from mint.devices_mi import *
import logging
ilename="logs/main.log"
logging.basicConfig( level=logging.DEBUG)
logger = logging.getLogger(__name__)


tws_i1 = Twiss()
tws_i1.E = 0.005000000
tws_i1.beta_x  = 55.7887190242
tws_i1.beta_y  = 55.7887190242
tws_i1.alpha_x = 18.185436973
tws_i1.alpha_y = 18.185436973
tws_i1.beta_x  = 0.286527307369
tws_i1.beta_y  = 0.286527307369
tws_i1.alpha_x = -0.838833736086
tws_i1.alpha_y = -0.838833736086

class BigCavity:
    def __init__(self, eid):
        self.id = eid
        self.cavs = []

    @property
    def phi(self):
        return np.mean([cav.phi for cav in self.cavs])

    @phi.setter
    def phi(self, value):
        for cav in self.cavs:
            cav.phi = value

    @property
    def v(self):
        return np.sum([cav.v for cav in self.cavs])

    @v.setter
    def v(self, value):
        cav_v = value/len(self.cavs)
        for cav in self.cavs:
            cav.v = cav_v

class TrackTread(Thread):
    def __init__(self, sections, phys_procs, check_header, gui):
        super(TrackTread, self).__init__()
        self.sections = sections
        self.phys_procs = phys_procs
        self.check_header = check_header
        self.gui = gui
        self.particles = None
        self.kill_track = False

    def kill_tracking(self):
        for sec in self.sections:
            sec.kill_track = True
            if "navigator" in sec.__dict__.keys():
                sec.navigator.kill_process = True

    def start_tracking(self):
        for sec in self.sections:
            sec.kill_track = False
            if "navigator" in sec.__dict__.keys():
                sec.navigator.kill_process = False

    def run(self):
        particles = self.particles
        t0 = time.time()
        self.start_tracking()
        for sec in self.sections:
            if sec.kill_track:
                continue
            for proc in self.phys_procs:
                flag = sec.ui.get_status(proc=proc.__name__)
                flag_name = sec.translator[proc]
                sec.__dict__[flag_name] = flag
                #print(proc.__name__, flag_name, flag)
            if sec.ui.get_status(proc=self.check_header):
                #print("tracking trough "+sec.__class__.__name__ + " ....")
                #sec.calc_tws = False
                particles = sec.tracking(particles)
                #sec.calculated = True
                sec.ui.calculated(True)
                self.gui.plot_s2e(particles)
            print(sec.__class__.__name__ + "time exec: ", time.time() - t0)


class S2ETool:
    """ Main class for the GUI application """
    def __init__(self):
        self.gui = OcelotInterfaceWindow(master=self)
        self.ui = self.gui.ui
        self.tws0 = tws_i1
        self.online_calc = True
        self.phys_procs = [CSR, SpaceCharge, WakeKick, BeamTransform]
        self.check_header = "Calculate?"
        self.load_sections()
        self.load_s2e_table()
        self.load_quads()
        cav_names = ["A1", "AH1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10",
                     "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19", "A20",
                     "A21", "A22", "A23", "A24", "A25"]
        self.create_super_cavs(cav_names)

        self.load_super_cavs()
        self.set_cavs()

        self.mi = XFELMachineInterface()
        self.magnets = MIMagnets()
        self.magnets.mi = self.mi
        self.mi_cav_obj = MICavity()
        self.mi_cav_obj.mi = self.mi
        self.track_thread = None

    def stop_s2e(self):
        if self.track_thread != None:
            self.track_thread.stop_flag = True
            self.track_thread.kill_tracking()

    def start_s2e(self):
        self.update_lattice_from_tables()
        particles = None
        for sec in self.sections:
            sec.ui.calculated(False)
        self.track_thread = TrackTread(self.sections, self.phys_procs, self.check_header, self.gui)
        self.track_thread.start()

    def set_cavs(self):
        for cav in self.cavs:
            if cav.id == "A1":
                v_i1 = 0.02118205
                phi_i1 = +28.4961 - 0.13965
                cav.v = v_i1*8
                cav.phi = phi_i1
            elif cav.id == "AH1":
                v3_i1 = 0.003552325 + 0.00004
                phi3_i1 = -153.5425 + 0.1 + 0 * 1.059098
                cav.v = v3_i1*8
                cav.phi = phi3_i1
            elif cav.id == "A2":
                v_l1 = 0.02005201875
                phi_l1 = 27.5615 - 0.0041
                cav.v = v_l1*32
                cav.phi = phi_l1
            elif cav.id in ["A3", "A4", "A5"]:
                v_l2 = 0.01913229167
                phi_l2 = 21.8576 + 0.0833
                cav.v = v_l2*32
                cav.phi = phi_l2
            elif cav.id in ["A6", "A7", "A8", "A9", "A10",
                                      "A11", "A12", "A13", "A14", "A15", "A16", "A17","A18", "A19", "A20",
                                      "A21", "A22", "A23", "A24", "A25"]:
                v_l3 = 0.02248065476
                phi_l3 = 0.0
                cav.v = v_l3*32
                cav.phi = phi_l3
            cav.ui.set_volt(cav.v)
            cav.ui.set_phi(cav.phi)


    def load_s2e_table(self):
        all_phys_procs = [proc.__name__ for proc in self.phys_procs]

        self.gui.init_s2e_table(self.sections, self.check_header, phys_proc_names=all_phys_procs)

    def load_quads(self):
        self.quads = []
        for sec in self.sections:
            for elem in sec.lattice.sequence:
                if elem.__class__ == Quadrupole:
                    self.quads.append(elem)

        self.gui.init_quad_table(self.quads, calc_obj=self.calc_twiss)

    def update_lat(self):
        for sec in self.sections:
            sec.lattice.update_transfer_maps()

    def create_super_cavs(self, cav_names):
        self.cavs = []
        for cav_name in cav_names:
            cav = BigCavity(eid=cav_name)
            self.cavs.append(cav)


    def load_super_cavs(self):
        for sec in self.sections:

            for elem in sec.lattice.sequence:
                if elem.__class__ == Cavity:
                    for cav in self.cavs:
                        if "." + cav.id + "." in elem.id:
                            cav.cavs.append(elem)

                    #self.cavs.append(elem)

        self.gui.init_cavs_table(self.cavs, calc_obj=self.calc_twiss)

    def calc_twiss(self, calc=True):
        #lat = MagneticLattice(cell)
        if self.online_calc == False:
            return

        # L = 0
        tws0 = deepcopy(self.tws0)
        tws_all = []
        self.table2cavs()
        self.table2quads()
        for sec in self.sections:
        #    for elem in sec.lattice.sequence:
        #        if elem.__class__ in [Quadrupole]:
        #            #print(elem.id, elem.row)
        #            #elem.kick_mrad = elem.ui.get_value()
        #            elem.k1 = elem.ui.get_value()
        #            elem.ui.value_was_changed(np.abs(np.abs(elem.ui.get_init_value()) - np.abs(elem.k1))> 0.1)
        #
        #
        #                #self.r_items[elem.ui.row].setBrush(pg.mkBrush("r"))
        #                #self.ui.tableWidget.item(elem.row, 1).setForeground(QtGui.QColor(255, 101, 101))  # red
        #            #else:
        #                #self.ui.tableWidget.item(elem.row, 1).setForeground(QtGui.QColor(255, 255, 255))  # white
        #                #self.r_items[elem.ui.row].setBrush(pg.mkBrush("g"))
        #            #r = self.r_items[elem.ui.row]
        #            #sizes = r.init_params
        #            #sizes = list(r.boundingRect().getRect())
        #            #sizes[3] = 10*elem.kick_mrad/self.quad_ampl
        #            #r.setRect(sizes[0], sizes[1], sizes[2], sizes[3])
        #
            sec.lattice.update_transfer_maps()
            tws = twiss(sec.lattice, tws0)
            tws0 = tws[-1]
            tws_all = np.append(tws_all, tws)

        beta_x = [tw.beta_x for tw in tws_all]
        beta_y = [tw.beta_y for tw in tws_all]
        dx = [tw.Dx for tw in tws_all]
        dy = [tw.Dy for tw in tws_all]
        s = [tw.s for tw in tws_all]

        self.ui.w_twiss_monitor.update_plot(s, beta_x, beta_y, dx, dy)

    def table2quads(self):
        for elem in self.quads:
            k1 = elem.ui.get_value()
            elem.k1 = k1
            elem.ui.value_was_changed(np.abs(np.abs(elem.ui.get_init_value() - k1)) > 0.1)

    def table2cavs(self):
        for elem in self.cavs:
            v = elem.ui.get_volt()
            phi = elem.ui.get_phi()
            elem.v = v
            elem.phi = phi
            elem.ui.phi_was_changed(np.abs(np.abs(elem.ui.get_init_phi() - phi)) > 0.01)
            elem.ui.volt_was_changed(np.abs(np.abs(elem.ui.get_init_volt() - v)) > 0.01)
            #print(elem.id , v, phi)
            #elem.transfer_map = self.lat.method.create_tm(elem)

    def update_lattice_from_tables(self):
        self.table2quads()
        self.table2cavs()
        self.update_lat()

    def read_from_doocs(self):
        self.online_calc = False
        self.magnets.get_magnets(self.quads)
        self.mi_cav_obj.get_cavities(self.cavs)
        self.update_lattice_from_tables()
        self.online_calc = True
        self.calc_twiss()
        print("Calculate")

    def init_lattice(self):
        pass

    #def calc_twiss(self):
    #    pass

    def part_tracking(self):
        pass

    def load_particles(self):
        pass

    def estimate_fel(self):
        pass

    def load_sections(self):
        blank_sections = sections.sections
        self.sections = []
        for sec in blank_sections:
            self.sections.append(sec())
        #for sec in self.sections:
        #    print(sec.__class__.__name__)
        return self.sections

    def get_whole_lat(self):
        pass





def main():

    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)

    #create the application
    app = QApplication(sys.argv)
    # setting the path variable for icon
    #path = os.path.join(os.path.dirname(sys.modules[__name__].__file__), 'ocelot.png')
    #app.setWindowIcon(QtGui.QIcon(path))
    window = S2ETool()


    #window.setWindowIcon(QtGui.QIcon('ocelot.png'))
    window.gui.show()

    #Build documentaiton if source files have changed
    # TODO: make more universal
    #os.system("cd ./docs && xterm -T 'Ocelot Doc Builder' -e 'bash checkDocBuild.sh' &")
    #exit script
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
