"""
Online optics monitor.
S.Tomin. XFEL/DESY. 2017
"""
import sys, os
path = os.path.realpath(__file__)
indx = path.find("gui_tools")
print("PATH", os.path.realpath(__file__), path[:indx])
sys.path.append(path[:indx])
from ocelot import *
from ocelot.cpbd.magnetic_lattice import shrinker
#from accelerator.s2e_sections import sections
from accelerator.lattice.xfel_i1_mad import *
from accelerator.lattice.xfel_l1_mad import *
from accelerator.lattice.xfel_l2_mad import *
from accelerator.lattice.xfel_l3_no_cl_mode_B import *
from accelerator.lattice.xfel_cl_mode_B import *
from accelerator.lattice.xfel_sase1_mode_B import *
from accelerator.lattice.xfel_sase3_mode_B import *
from mint.devices_mi import *
from gui_tools.gui.gui_optics_mon import OcelotOpticsWindow
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from PyQt5 import QtCore
from PyQt5.QtWidgets import QMainWindow, QApplication
from ocelot.optimizer.mint.xfel_interface import *
import numpy as np
from copy import deepcopy
tws_i1 = Twiss()
tws_i1.E = 0.005000000
tws_i1.beta_x  = 55.7887190242
tws_i1.beta_y  = 55.7887190242
tws_i1.alpha_x = 18.185436973
tws_i1.alpha_y = 18.185436973

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

class OpticsMonitor:
    """ Main class for the GUI application """
    def __init__(self):
        self.gui = OcelotOpticsWindow(master=self)
        self.ui = self.gui.ui
        self.tws0 = tws_i1
        #self.load_sections()
        #self.load_s2e_table()
        self.define_lattice()
        self.load_quads()
        #self.load_cavities(cav_names=["A1", "AH1", "A2", "A3"])
        self.mi = XFELMachineInterface()
        self.magnets = MIMagnets()
        self.magnets.mi = self.mi
        self.online_calc = True
        self.tws_meas = Twiss()
        self.tws0 = tws_i1

    def auto_reading(self):
        self.online_calc = False
        self.magnets.get_magnets(self.quads)
        self.table2quads()
        self.online_calc = True
        self.back_tracking()
        self.calc_twiss()
        print("Calculate")

    def table2quads(self):
        for elem in self.quads:
            k1 = elem.ui.get_value()
            elem.k1 = k1
            elem.transfer_map = self.lat.method.create_tm(elem)

    def load_quads(self):
        self.quads = []
        for elem in self.lat.sequence:
            if elem.__class__ == Quadrupole:
                self.quads.append(elem)
        self.gui.init_quad_table(self.quads, calc_obj=self.calc_twiss)

    def define_lattice(self):
        self.cell_back_track = (cell_i1 + cell_l1 + cell_l2 + cell_l3_no_cl + cell_cl)
        seq = cell_i1 + cell_l1 + cell_l2 + cell_l3_no_cl + cell_cl + cell_sase1 + cell_t4 + cell_sase3
        self.lat = MagneticLattice(seq)
        self.lat = shrinker(self.lat, remaining_types=[Quadrupole, Bend, SBend, RBend, Cavity], 
        remaining_elems = ["OTRC.55.I1", "OTRB.218.B1", "OTRB.450.B2"],
                            init_energy=self.tws0.E)

    def calc_twiss(self):
        if self.online_calc == False:
            return
        #tws0 = deepcopy(self.tws0)

        for elem in self.lat.sequence:
            if elem.__class__ in [Quadrupole]:
                #print(elem.id, elem.row)
                #elem.kick_mrad = elem.ui.get_value()
                elem.k1 = elem.ui.get_value()
                elem.ui.value_was_changed(np.abs(np.abs(elem.ui.get_init_value()) - np.abs(elem.k1))> 0.1)

                elem.transfer_map = self.lat.method.create_tm(elem)
        #print(self.tws0)
        tws = twiss(self.lat, self.tws0)


        beta_x = [tw.beta_x for tw in tws]
        beta_y = [tw.beta_y for tw in tws]
        dx = [tw.Dx for tw in tws]
        dy = [tw.Dy for tw in tws]
        s = [tw.s for tw in tws]

        self.ui.w_twiss_monitor.update_plot(s, beta_x, beta_y, dx, dy)

    def read_twiss(self):
        tws = Twiss()
        if self.ui.cb_otrb218.isChecked():
            section = "B1"
            tws.E = 0.7
        elif self.ui.cb_otrb450.isChecked():
            section = "B2"
            tws.E = 2.4
        elif self.ui.cb_otrc55.isChecked():
            #self.ui.cb_otrc55.setChecked(True)
            section = "I1"
            tws.E = 0.130
        else:
            tws.beta_x = 2.36238404123
            tws.beta_y = 2.90712039319
            tws.alpha_x = 1.23079453323
            tws.alpha_y = -1.45354874634
            tws.E = 0.130
            return tws

        ch_beta_x = "XFEL.UTIL/BEAM_PARAMETER/" + section + "/PROJECTED_X.BETA.SA1"
        ch_alpha_x = "XFEL.UTIL/BEAM_PARAMETER/" + section + "/PROJECTED_X.ALPHA.SA1"
        ch_beta_y = "XFEL.UTIL/BEAM_PARAMETER/" + section + "/PROJECTED_Y.BETA.SA1"
        ch_alpha_y = "XFEL.UTIL/BEAM_PARAMETER/" + section + "/PROJECTED_Y.ALPHA.SA1"
        ch_energy = "XFEL.UTIL/BEAM_PARAMETER/" + section + "/PROJECTED_X.ENERGY.SA1"

        tws.beta_x = self.mi.get_value(ch_beta_x)
        tws.beta_y = self.mi.get_value(ch_beta_y)
        tws.alpha_x = self.mi.get_value(ch_alpha_x)
        tws.alpha_y = self.mi.get_value(ch_alpha_y)
        #tws.E = self.mi.get_value(ch_energy)*0.001
        print(tws)
        return tws

    def back_tracking(self):
        tws0 = self.read_twiss()
        if np.abs(tws0.beta_x - self.tws_meas.beta_x)>0.01 or np.abs(tws0.beta_y - self.tws_meas.beta_y)>0.01:
            self.tws_meas = tws0
        else:
            return self.tws0
        #if self.ui.cb_design_tws.isChecked():
        #    return self.tws_des
        if self.ui.cb_otrb218.isChecked():
            stop = otrb_218_b1

        elif self.ui.cb_otrb450.isChecked():
            stop = otrb_450_b2

        else:
            stop = otrc_55_i1


        lat_tmp = MagneticLattice(self.lat.sequence, stop=stop)
        lat_tmp = MagneticLattice(lat_tmp.sequence[::-1])
        for elem in lat_tmp.sequence:
            if elem.__class__ == Quadrupole:
                print(elem.id, elem.k1)
            if elem.__class__  == Cavity:
                elem.phi -= 180
                #print(elem.v, elem.phi)
                elem.transfer_map = self.lat.method.create_tm(elem)
        #lat_tmp.update_transfer_maps()

        tws0.alpha_x = -tws0.alpha_x
        tws0.alpha_y = -tws0.alpha_y
        #print("start", tws0)
        tws = twiss(lat_tmp, tws0)
        #plot_opt_func(lat_tmp, tws)
        #plt.show()
        self.tws0 = Twiss()
        self.tws0.beta_x = tws[-1].beta_x
        self.tws0.beta_y = tws[-1].beta_y
        self.tws0.alpha_x = -tws[-1].alpha_x
        self.tws0.alpha_y = -tws[-1].alpha_y
        self.tws0.s = 0
        self.tws0.E = tws[-1].E
        for elem in lat_tmp.sequence:
            if elem.__class__  == Cavity:
                elem.phi += 180
                elem.transfer_map = self.lat.method.create_tm(elem)
        lat_tmp.update_transfer_maps()
        return self.tws0

def main():

    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)

    #create the application
    app = QApplication(sys.argv)
    # setting the path variable for icon
    #path = os.path.join(os.path.dirname(sys.modules[__name__].__file__), 'ocelot.png')
    #app.setWindowIcon(QtGui.QIcon(path))
    window = OpticsMonitor()


    #window.setWindowIcon(QtGui.QIcon('ocelot.png'))
    window.gui.show()

    #Build documentaiton if source files have changed
    # TODO: make more universal
    #os.system("cd ./docs && xterm -T 'Ocelot Doc Builder' -e 'bash checkDocBuild.sh' &")
    #exit script
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
