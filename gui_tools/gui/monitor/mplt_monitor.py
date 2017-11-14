"""
S.Tomin Matplotlib widget
"""

import sys, os
from os import listdir
from os.path import isfile, join

from PyQt5.QtWidgets import QCheckBox,QSizePolicy, QPushButton, QHBoxLayout, QHeaderView, QApplication, QMenu, QWidget, QAction, \
    QTableWidget, QTableWidgetItem, QDoubleSpinBox, QGridLayout

from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from gui.monitor.ui_prtcl_monitor import Ui_Widget
import numpy as np
import pyqtgraph as pg
#from gui.monitor.ocl_monitor import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from ocelot.cpbd.beam import global_slice_analysis_extended
from ocelot import *
from scripts.io import *


class MpltMonitor(QWidget):
    def __init__(self, parent=None):
        super().__init__()
        self.ui = Ui_Widget()
        self.ui.setupUi(self)
        self.ui.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.ui.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.add_plot()
        self.loadStyleSheet()
        self.ui.pb_plot.clicked.connect(self.plot)
        #mypath = "C:/Users/tomins/Documents/Dropbox/DESY/repository/ocelot/gui_tools/track_data/particles"
        path = os.path.realpath(__file__)
        indx = path.find("gui_tools")
        self.p_dict = path[:indx] + "gui_tools"+ os.sep + "track_data"+os.sep + "particles"+os.sep
        self.add_files(self.p_dict)
        #print("PATH",self.p_dict) # os.path.realpath(__file__), path[:indx])
        self.p_arrays = {}
        self.slices = {}
        self.p_array = None
        self.slice_params = None
        #sys.path.append(path[:indx])
        self.ui.cb_p_distrib.addItem("E-S")
        self.ui.cb_p_distrib.addItem("Y-X")
        self.ui.cb_p_distrib.addItem("Y-S")
        self.ui.cb_p_distrib.addItem("X-S")
        self.ui.cb_p_distrib.addItem("Px-X")
        self.ui.cb_p_distrib.addItem("Py-Y")
        #self.ui.cb_p_distrib.setCurrentIndex(0)
        self.ui.cb_p_distrib.currentIndexChanged.connect(self.update_distrib_plot)

        self.ui.cb_slice_params.addItem("Current")
        self.ui.cb_slice_params.addItem("Emittance")
        self.ui.cb_slice_params.addItem("Energy")
        self.ui.cb_slice_params.addItem("Energy Spread")
        #self.ui.cb_slice_params.addItem("Px-X")
        #self.ui.cb_slice_params.addItem("Py-Y")
        self.ui.cb_slice_params.currentIndexChanged.connect(self.update_slice_plot)
        self.ui.pb_plot.clicked.connect(self.plot)
        self.ui.pb_reload_files.clicked.connect(self.reload_files)

    def reload_files(self):
        self.add_files(self.p_dict)
        self.p_arrays = {}
        self.slices = {}
        self.p_array = None
        self.slice_params = None

    def add_files(self, p_path):
        self.p_files = [f for f in listdir(p_path) if isfile(join(p_path, f))]
        #self.ui.cb_p_file.currentIndexChanged.connect(lambda: 0)
        try:
            self.ui.cb_p_file.currentIndexChanged.disconnect()
        except:
            pass
        self.ui.cb_p_file.clear()
        for path in self.p_files:
            self.ui.cb_p_file.addItem(path)
        #self.ui.cb_p_file.setCurrentIndex(0)
        self.ui.cb_p_file.currentIndexChanged.connect(self.plot)

    def calculate_slaice_params(self, p_array):
        slice_params = global_slice_analysis_extended(p_array, 5000, 0.01, 2, 2)
        return slice_params

    def load_particles(self):
        current_file = self.ui.cb_p_file.currentText()
        #if current_file in self.p_arrays.keys() and current_file in self.slices.keys():
        #    print("Already loaded")
        #    return
        self.p_arrays[current_file] = read_beam_file(self.p_dict + current_file)
        self.slices[current_file] = self.calculate_slaice_params(self.p_arrays[current_file])
        #s, I, ex, ey, me, se, gamma0, emitxn, emityn = slice_params

    def update_slice_plot(self):
        #s, I, ex, ey, me, se, gamma0, emitxn, emityn = self.slice_params
        #print(self.slice_params[6:])
        current_slice_param = self.ui.cb_slice_params.currentText()
        self.ax_l.clear()
        if current_slice_param == "Current":
            self.ax_l.plot(self.slice_params[0]*1e3, self.slice_params[1], "r", label="Current")
            self.ax_l.set_xlabel("s, mm")
            self.ax_l.set_ylabel("I, A")
        elif current_slice_param == "Emittance":
            self.ax_l.plot(self.slice_params[0]*1e3, self.slice_params[2], "r", label="emit_x")
            self.ax_l.plot(self.slice_params[0]*1e3, self.slice_params[3], "b", label="emit_y")
            self.ax_l.set_xlabel("s, mm")
            self.ax_l.set_ylabel("emit, mm*mrad")
        elif current_slice_param == "Energy":
            self.ax_l.plot(self.slice_params[0]*1e3, self.slice_params[4], "r", label="Energy")
            self.ax_l.set_xlabel("s, mm")
            self.ax_l.set_ylabel("E, eV")
        else:
            self.ax_l.plot(self.slice_params[0]*1e3, self.slice_params[5], "b", label="Energy Spread")
            self.ax_l.set_xlabel("s, mm")
            self.ax_l.set_ylabel("dE, eV")
        self.ax_l.grid(True)
        self.ax_l.legend()
        self.canvas.draw()


    def update_distrib_plot(self):
        current_dist = self.ui.cb_p_distrib.currentText()
        #print(self.p_array.x())
        self.ax_r.clear()
        if current_dist == "E-S":
            self.ax_r.plot(self.p_array.tau()*1e3, self.p_array.p(), "r.", label=current_dist)
            self.ax_r.set_xlabel("S, mm")
            self.ax_r.set_ylabel("dE/E")
        elif current_dist == "Y-X":
            self.ax_r.plot(self.p_array.x()*1e3, self.p_array.y()*1e3, "b.", label=current_dist)
            self.ax_r.set_xlabel("X, mm")
            self.ax_r.set_ylabel("Y, mm")
        elif current_dist == "Y-S":
            self.ax_r.plot(self.p_array.tau()*1e3, self.p_array.y()*1e3, "b.", label=current_dist)
            self.ax_r.set_xlabel("tau, mm")
            self.ax_r.set_ylabel("Y, mm")
        elif current_dist == "X-S":
            self.ax_r.plot(self.p_array.tau()*1e3, self.p_array.x()*1e3, "b.", label=current_dist)
            self.ax_r.set_xlabel("tau, mm")
            self.ax_r.set_ylabel("Y, mm")
        elif current_dist == "Px-X":
            self.ax_r.plot(self.p_array.px(), self.p_array.x()*1e3, "b.", label=current_dist)
            self.ax_r.set_xlabel("px/p0")
            self.ax_r.set_ylabel("X, mm")
        else:
            self.ax_r.plot(self.p_array.py(), self.p_array.y()*1e3, "b.", label=current_dist)
            self.ax_r.set_xlabel("py/p0")
            self.ax_r.set_ylabel("Y, mm")
        self.ax_r.grid(True)
        self.ax_r.legend()
        self.canvas.draw()

    def add_plot(self):
        #plt.style.use('dark_background')
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        #self.button = QPushButton('Plot')
        #self.button.clicked.connect(self.plot)

        # set the layout
        layout = QGridLayout()

        self.ui.w_monitor.setLayout(layout)
        layout.addWidget(self.toolbar, 0, 0)
        layout.addWidget(self.canvas, 1, 0)
        #layout.addWidget(self.button, 2, 0)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

        self.ax_l = self.figure.add_subplot(121)
        self.ax_r = self.figure.add_subplot(122)
        self.ax_l.grid(True)
        self.ax_r.grid(True)
        #self.figure.tight_layout()

    def plot(self):
        current_file = self.ui.cb_p_file.currentText()
        if current_file not in self.p_arrays.keys() or current_file not in self.slices.keys():
            self.load_particles()

        self.p_array = self.p_arrays[current_file]
        self.slice_params = self.slices[current_file]
        self.update_slice_plot()
        self.update_distrib_plot()


    def plot_s2e(self, particles):
        # create an axis
        self.p_array = particles
        self.slice_params = self.calculate_slaice_params(particles)
        self.update_slice_plot()
        self.update_distrib_plot()


    def loadStyleSheet(self):
        """ Load in the dark theme style sheet. """
        try:
            self.cssfile = "style.css"
            with open(self.cssfile, "r") as f:
                self.setStyleSheet(f.read())
        except IOError:
            print('No style sheet found!')

        # @pyqtSlot()
        # def on_click(self):
        #    print("\n")
        #    for currentQTableWidgetItem in self.ui.tableWidget.selectedItems():
        #        print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MpltMonitor()
    window.show()
    sys.exit(app.exec_())