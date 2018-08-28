"""
This class takes care about GUI functionality.
S.Tomin, 2018
"""

from gui.UIOnline import *
#from table_test import *
import sys
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QAction, QTableWidget, QTableWidgetItem, QVBoxLayout, QFileDialog
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from PyQt5 import QtCore
from gui.device.device_ui import *
from ocelot import *
import os
import importlib.util
import numpy as np

class OcelotInterfaceWindow(QMainWindow):
    """ Main class for the GUI application """
    def __init__(self, master=None):
        super().__init__()

        # initialize
        #QFrame.__init__(self)

        #self.logbook = "xfellog"
        #self.dev_mode = False
        self.lattice_dir_root = "accelerator/lattice"
        self.master = master
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.menuBar.setNativeMenuBar(False)
        self.ui.mainToolBar.setVisible(False)
        self.loadStyleSheet()
        self.ui.gridLayout_6.setContentsMargins(0,0,0,0)

        self.ui.w_q_table.ui.pb_calculate.clicked.connect(self.master.calc_twiss)

        self.ui.pb_start_tracking.clicked.connect(self.ui_start_s2e)
        self.ui.pb_hide_s2e_pan.clicked.connect(self.hide_show_s2e)
        self.ui.pb_hide_cav_pan.clicked.connect(self.hide_show_cav)
        self.ui.horizontalLayout.setContentsMargins(0,0,0,0)
        self.ui.gridLayout_3.setContentsMargins(0,0,0,0)

        self.ui.pb_twiss_from_track.clicked.connect(self.plot_track_twiss)

        self.ui.pb_hide_show_quads.clicked.connect(self.hide_show_quads)
        self.ui.pb_read.clicked.connect(self.master.read_from_doocs)
        #self.stop_track_timer = QtCore.QTimer()
        #self.stop_track_timer.timeout.connect(self.check_tread_is_alive)

        self.ui.actionLoad_Lattice.triggered.connect(self.load_lattice)
        self.ui.actionSave_Lattice.triggered.connect(self.save_lattice)


    def load_lattice(self):
        filename = QFileDialog.getExistingDirectory(self, 'Load Lattice',
                                               self.lattice_dir_root,
                                               QFileDialog.DontUseNativeDialog)

        if filename:
            print(filename)
            for sec in self.master.sections:
                new_cell = {}
                fname = filename + os.sep + sec.lattice_name + ".py"


                spec = importlib.util.spec_from_file_location(sec.lattice_name, fname)
                foo = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(foo)
                #print(foo.cell)
                for elem in foo.cell:

                    new_cell[elem.id] = elem
                #print(self.master.quads)
                for elem in self.master.quads:
                    #print(elem.id, elem.id in new_cell.keys())
                    if elem.id in new_cell.keys():
                        elem.ui.set_value(new_cell[elem.id].k1)
                for elem in self.master.bends:
                    if elem.id in new_cell.keys():
                        #print("bend", elem.id, new_cell[elem.id].angle)
                        elem.ui.set_value(new_cell[elem.id].angle)
                for elem in self.master.cavs:
                    #print(elem.cavs[0].id, elem.cavs[0].id in new_cell.keys())
                    if elem.cavs[0].id in new_cell.keys():
                        elem.ui.set_volt(new_cell[elem.cavs[0].id].v*len(elem.cavs))
                        elem.ui.set_phi(new_cell[elem.cavs[0].id].phi)

                        #print("set cav: ", elem.id, new_cell[elem.cavs[0].id].v*len(elem.cavs))
                        #print("set cav: ", elem.id, new_cell[elem.cavs[0].id].phi)
            self.master.update_lattice_from_tables()

    def save_lattice(self):
        #print(self.master.gold_orbits_dir)
        filename = QFileDialog.getExistingDirectory(self, 'Save Lattice',
                                               self.lattice_dir_root,
                                               QFileDialog.DontUseNativeDialog)
        print(filename)

        if filename:
            for sec in self.master.sections:
                fname = filename + os.sep + sec.lattice_name + ".py"
                write_lattice(sec.lattice, file_name=fname, remove_rep_drifts=True, power_supply=False)

        #    name = filename.split("/")[-1]
        #    parts = name.split(".")
        #    body_name = parts[0]
        #
        #    if len(parts) < 2 or parts[1] != "json":
        #        part = filename.split(".")[0]
        #        filename = part + ".json"
        #
        #    self.save_golden_orbit(filename)

    def ui_start_s2e(self):
        if self.ui.pb_start_tracking.text() == "Stop Tracking":
            self.ui.pb_start_tracking.setText("Start Tracking")
            self.ui.pb_start_tracking.setStyleSheet("color: rgb(85, 255, 127);")
            self.master.stop_s2e()
        else:
            self.master.start_s2e()
            self.ui.pb_start_tracking.setText("Stop Tracking")
            self.ui.pb_start_tracking.setStyleSheet("color: rgb(255, 0, 0);")

    def closeEvent(self, event):
        self.ui.w_s2e.close()
        QMainWindow.closeEvent(self, event)


    def hide_show_quads(self):
        if self.ui.pb_hide_show_quads.text() == "Hide Quads Panel":
            self.ui.pb_hide_show_quads.setText("Show Quads Panel")
            self.ui.pb_hide_show_quads.setStyleSheet("color: rgb(255, 0, 255);")
            self.ui.w_q_table.hide()
        else:
            self.ui.w_q_table.show()
            self.ui.pb_hide_show_quads.setText("Hide Quads Panel")
            self.ui.pb_hide_show_quads.setStyleSheet("color: rgb(85, 255, 255);")

    def hide_show_s2e(self):
        if self.ui.pb_hide_s2e_pan.text() == "Hide S2E Panel":
            self.ui.pb_hide_s2e_pan.setText("Show S2E Panel")
            self.ui.pb_hide_s2e_pan.setStyleSheet("color: rgb(255, 0, 255);")
            self.ui.w_s2e.hide()
        else:
            self.ui.w_s2e.show()
            self.ui.pb_hide_s2e_pan.setText("Hide S2E Panel")
            self.ui.pb_hide_s2e_pan.setStyleSheet("color: rgb(85, 255, 255);")

    def hide_show_cav(self):
        if self.ui.pb_hide_cav_pan.text() == "Hide Cavity Panel":
            self.ui.pb_hide_cav_pan.setText("Show Cavity Panel")
            self.ui.pb_hide_cav_pan.setStyleSheet("color: rgb(255, 0, 255);")
            self.ui.w_cav_table.hide()
        else:
            self.ui.w_cav_table.show()
            self.ui.pb_hide_cav_pan.setText("Hide Cavity Panel")
            self.ui.pb_hide_cav_pan.setStyleSheet("color: rgb(85, 255, 255);")


    #def set_s2e_calc_obj(self, calc_obj):
    #    self.s2e_calc_obj = calc_obj
    def load_all_twiss(self):
        bx = []
        by = []
        s = []
        E = []
        for sec in self.master.sections:
            print(sec.lattice_name, os.path.isfile(sec.tws_file) )#if sec.tws_file
            if os.path.isfile(sec.tws_file):
                tws_data = sec.load_twiss_file()
            else:
                print(sec.lattice_name, "NO TWISS FILE")

            bx = np.append(bx, tws_data["beta_x"])
            by = np.append(by, tws_data["beta_y"])
            s0 = 0 if len(s) == 0 else s[-1]
            s = np.append(s, tws_data["s"]+s0)
            E = np.append(E, tws_data["E"])
        return s, bx, by, E


    def plot_track_twiss(self):
        s, bx, by, E = self.load_all_twiss()
        self.ui.w_twiss_monitor.update_plot_track(s, bx, by, E)

    def plot_s2e(self, particles):
        self.ui.w_s2e_monitor.plot_s2e(particles)

    def init_quad_table(self, devs, calc_obj):
        self.ui.w_q_table.init_table(devs=devs, device_iface=DeviceUI, calc_obj=calc_obj, spin_params=[-5000, 5000, 0.1])

    def init_bend_table(self, devs, calc_obj):
        for row, dev in enumerate(devs):
            ui = VirtualUI()
            #ui.tableWidget = self.ui.tableWidget
            ui.row = row
            ui.col = 2
            ui.set_value(dev.angle)
            ui.set_init_value(dev.angle)

            dev.ui = ui

    def init_s2e_table(self, sections, check_header="", phys_proc_names=[]):
        self.ui.w_s2e.init_table(sections=sections, check_header=check_header,
                                 phys_proc_names=phys_proc_names)

    def init_cavs_table(self, devs, calc_obj):
        self.ui.w_cav_table.init_table(devs=devs, device_iface=CavityUI, calc_obj=calc_obj,
                                       spin_params_phi=[-360, 360, 0.1],
                                       spin_params_v=[0, 1000, 0.001])


    def loadStyleSheet(self):
        """ Load in the dark theme style sheet. """
        try:
            self.cssfile = "gui/style.css"
            with open(self.cssfile, "r") as f:
                self.setStyleSheet(f.read())
        except IOError:
            print('No style sheet found!')

def main():

    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)

    #create the application
    app = QApplication(sys.argv)
    # setting the path variable for icon
    #path = os.path.join(os.path.dirname(sys.modules[__name__].__file__), 'ocelot.png')
    #app.setWindowIcon(QtGui.QIcon(path))
    window = OcelotInterfaceWindow()


    #window.setWindowIcon(QtGui.QIcon('ocelot.png'))
    window.show()

    #Build documentaiton if source files have changed
    # TODO: make more universal
    #os.system("cd ./docs && xterm -T 'Ocelot Doc Builder' -e 'bash checkDocBuild.sh' &")
    #exit script
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
