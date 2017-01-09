#QT imports
from __future__ import absolute_import, print_function

import sys
import os
path = sys.path[0]
indx = path.find("ocelot")
sys.path.append(path[:indx])
# for pyqtgraph import
sys.path.append(path[:indx]+"ocelot")

from PyQt4.QtGui import QApplication, QFrame, QPixmap, QMessageBox
from PyQt4 import QtGui, QtCore

#normal imports
import numpy as np
import subprocess
import platform
import time
import pyqtgraph as pg
if sys.version_info[0] == 2:
    from imp import reload
else:
    from importlib import reload

##logbook imports
#from re import sub
#from xml.etree import ElementTree
#from shutil import copy
#from datetime import datetime
#import Image
#import json
#GUI layout file
#from ocelot.optimizer.UIOcelotInterface_gen import Ui_Form
from ocelot.optimizer.gui_main import *
#slac python toolbox imports
#from ocelot.optimizer.mint.lcls_interface import MatLog
#from ocelot.optimizer.mint.lcls_interface import TestMatLog as MatLog
#local imports

#from ocelot.optimizer import scanner_threads
from ocelot.optimizer.mint.opt_objects import *
from ocelot.optimizer.mint import mint
from ocelot.optimizer.mint import opt_objects as obj
from ocelot.optimizer.mint import obj_function
from ocelot.utils import db
from ocelot.optimizer.mint.xfel_interface import *

class OcelotInterfaceWindow(QFrame):
    """ Main class for the GUI application """
    def __init__(self):
        """
        Initialize the GUI and QT UI aspects of the application.

        Create the epicsGet class, to try and get around network errors.
        Initialize the scan parameters.
        Connect start and logbook buttons on the scan panel.
        Initialize the plotting.
        Make the timer object that updates GUI on clock cycle durring a scan.
        """

        self.set_file = "./parameters/default.json"
        # initialize
        QFrame.__init__(self)

        # try
        self.ui = MainWindow(self)

        #self.mi = TestMachineInterface()
        self.mi = XFELMachineInterface()
        self.dp = TestDeviceProperties(ui=self.ui.widget)

        self.opt_control = mint.OptControl()
        self.objective_func = obj_function.XFELTarget()
        self.objective_func_pv = "test_obj"

        self.addPlots()

        # database
        self.dbname = "./parameters/test.db"
        # db.create_db(self.dbname)
        self.db = db.PerfDB(dbname=self.dbname)

        #object funciton selectinator (gdet)
        #self.setObFunc()

        #load in the dark theme style sheet
        #self.loadStyleSheet()
        self.hyper_file = "../parameters/hyperparameters.npy"

        self.ui.pb_start_scan.clicked.connect(self.start_scan)
        self.ui.pb_edit_obj_func.clicked.connect(self.run_editor)
        self.ui.cb_use_predef.stateChanged.connect(self.set_obj_fun)
        self.ui.pb_logbook.clicked.connect(self.ui.logbook)

        #self.opt = mint.Optimizer()
        self.ui.restore_state(self.set_file)
        self.path_to_obj_func = os.path.join(os.path.dirname(sys.modules[__name__].__file__), 'mint/obj_function.py')


        self.set_obj_fun()
        self.m_status = mint.MachineStatus()
        self.set_m_status()
        self.opt_control = mint.OptControl()
        self.opt_control.m_status = self.m_status

        self.name1 = "Nelder-Mead Simplex"
        self.name2 = "Gaussian Process"
        self.name3 = "Custom Minimizer"
        self.name4 = "Conjugate Gradient"
        self.name5 = "Powell's Method"
        self.ui.cb_select_alg.addItem(self.name1)
        self.ui.cb_select_alg.addItem(self.name2)
        self.ui.cb_select_alg.addItem(self.name3)
        #timer for plots, starts when scan starts
        self.multiPvTimer = QtCore.QTimer()
        self.multiPvTimer.timeout.connect(self.getPlotData)

        self.indicator = QtCore.QTimer()
        self.indicator.timeout.connect(self.indicate_machine_state)
        self.indicator.start(100)
        #self.ui.widget_2.setStyleSheet("background-color:red;")
        #self.ui.widget_3.setStyleSheet("background-color:red;")
        #p = w.palette()
        #p.setColor(w.backgroundRole(), QtCore.red)
        #w.setPalette(p)

    def scan_method_select(self):
        """
        Sets scanner method from options panel combo box selection.

        This method executes from the runScan() method, when the UI "Start Scan" button is pressed.

        Returns:
                 Selected scanner object
                 These objects are contrained in the scannerThreads.py file
        """
        index = self.ui.cb_select_alg.currentIndex()

        #GP Method
        if index == 1:
            minimizer = mint.GaussProcess()

        # Custom Minimizer
        elif index == 2:
            minimizer = mint.CustomMinimizer()

        # Conjugate Gradient
        #if index == 3:
        #    scanner = scanner_threads.OcelotScanner(parent=self, method='cg')
        #
        ## Powells Method
        #if index == 4:
        #    scanner = scanner_threads.OcelotScanner(parent=self,method='powell')

        #simplex Method
        else: #index == 0:
            minimizer = mint.Simplex()
        return minimizer


    def closeEvent(self, event):
        self.ui.save_state(self.set_file)
        if self.ui.pb_start_scan.text() == "Stop scan":
            self.opt.opt_ctrl.stop()
            del(self.opt)
            self.ui.pb_start_scan.setStyleSheet("color: rgb(85, 255, 127);")
            self.ui.pb_start_scan.setText("Start scan")
            return 0
        QFrame.closeEvent(self, event)

    def start_scan(self):

        self.scanStartTime = time.time()

        if self.ui.pb_start_scan.text() == "Stop scan":
            self.opt.opt_ctrl.stop()
            self.save2db()
            del(self.opt)
            self.ui.pb_start_scan.setStyleSheet("color: rgb(85, 255, 127);")
            self.ui.pb_start_scan.setText("Start scan")
            return 0

        self.pvs = self.ui.widget.getPvsFromCbState()
        self.devices = self.ui.widget.get_devices(self.pvs)

        if len(self.devices) == 0:
            self.error_box(message="Check Devices")
            return 0

        self.setUpMultiPlot(self.devices)
        self.multiPvTimer.start(100)

        self.set_obj_fun()
        #self.objective_func = obj_function.TestTarget()
        minimizer = self.scan_method_select()

        if minimizer.__class__ == mint.GaussProcess:
            minimizer.seed_iter = self.ui.sb_seed_iter.value()
            minimizer.seed_timeout = self.ui.sb_tdelay.value()
            minimizer.hyper_file = self.hyper_file

        elif minimizer.__class__ == mint.Simplex:
            if self.ui.cb_use_isim.checkState():
                minimizer.dev_steps = []

                for dev in self.devices:
                    if dev.simplex_step == 0:
                        lims = dev.get_limits()
                        rel_step = self.ui.sb_isim_rel_step.value()
                        minimizer.dev_steps.append((lims[1] - lims[0])*rel_step/100.)

            else:
                minimizer.dev_steps = None

        #minimizer = mint.Simplex()
        #minimizer = mint.CustomMinimizer()
        self.max_iter = self.ui.sb_num_iter.value()
        minimizer.max_iter = self.max_iter

        self.opt = mint.Optimizer()

        #self.opt_control = mint.OptControl()
        self.opt_control.m_status = self.m_status
        self.opt.opt_ctrl = self.opt_control
        self.opt.timeout = self.ui.sb_tdelay.value()

        self.opt.minimizer = minimizer

        seq = [mint.Action(func=self.opt.max_target_func, args=[self.objective_func, self.devices])]
        self.opt.seq = seq

        #self.opt.eval(seq)
        self.opt.start()
        self.ui.pb_start_scan.setText("Stop scan")
        self.ui.pb_start_scan.setStyleSheet("color: red")

    def scan_finished(self):
        try:
            if not self.opt.isAlive() and self.ui.pb_start_scan.text() == "Stop scan":
                self.ui.pb_start_scan.setStyleSheet("color: rgb(85, 255, 127);")
                self.ui.pb_start_scan.setText("Start scan")
                self.save2db()
        except:
            pass

    def indicate_machine_state(self):
        #print(self.opt_control.is_ok)
        if not self.opt_control.is_ok:
            self.ui.widget_2.setStyleSheet("background-color:red;")
            self.ui.widget_3.setStyleSheet("background-color:red;")
        else:    #time.sleep(0.5)
            self.ui.widget_2.setStyleSheet("background-color:323232;")
            self.ui.widget_3.setStyleSheet("background-color:323232;")

    def save2db(self):
        d_names = []
        d_start = []
        d_stop = []
        for dev in self.devices:
            d_names.append(dev.eid + "_val")
            d_start.append(dev.values[0])
            d_stop.append(dev.values[1])

            d_names.append(dev.eid + "_lim")
            d_start.append(dev.get_limits()[0])
            d_stop.append(dev.get_limits()[1])

        o_names = ["obj_id", "obj_value", "obj_pen", "niter"]
        o_start = [self.objective_func.eid, self.objective_func.values[0], self.objective_func.penalties[0], self.max_iter]
        o_stop = [self.objective_func.eid, self.objective_func.values[-1], self.objective_func.penalties[-1],  len(self.objective_func.penalties)]

        self.db.new_tuning({'wl': 13.6, 'charge': 0.1, 'comment': 'test tuning'})
        tune_id = self.db.current_tuning_id()

        start_sase = self.objective_func.values[0]
        stop_sase = self.objective_func.values[-1]
        self.db.new_action(tune_id, start_sase=start_sase, end_sase=stop_sase)

        #print('current actions in tuning', [(t.id, t.tuning_id, t.sase_start, t.sase_end) for t in self.db.get_actions()])

        action_id = self.db.current_action_id()

        self.db.add_action_parameters(tune_id, action_id, param_names=o_names+d_names, start_vals=o_start+d_start,
                                 end_vals=o_stop+d_stop)
        print('current actions', [(t.id, t.tuning_id, t.sase_start, t.sase_end) for t in self.db.get_actions()])
        print('current action parameters', [(p.tuning_id, p.action_id, p.par_name, p.start_value, p.end_value) for p in
                                            self.db.get_action_parameters(tune_id, action_id)])


    def is_le_addr_ok(self, line_edit):
        dev = str(line_edit.text())
        state = True
        try:
            self.mi.get_value(dev)
        except:
            state = False
        if state:
            line_edit.setStyleSheet("color: rgb(85, 255, 0);")
        else:
            line_edit.setStyleSheet("color: red")
        return state

    def set_obj_fun(self):

        self.ui.use_predef_fun()

        if self.ui.cb_use_predef.checkState():
            print("RELOAD")
            reload(obj_function)
        else:
            a_str = str(self.ui.le_a.text())
            self.is_le_addr_ok(self.ui.le_a)

            b_str = str(self.ui.le_b.text())
            self.is_le_addr_ok(self.ui.le_b)

            c_str = str(self.ui.le_c.text())
            self.is_le_addr_ok(self.ui.le_c)

            func = str(self.ui.le_obf.text())
            self.is_le_addr_ok(self.ui.le_obf)

            def get_value_exp():
                A = self.objective_func.mi.get_value(a_str)
                B = self.objective_func.mi.get_value(b_str)
                C = self.objective_func.mi.get_value(c_str)
                return eval(func)

        self.objective_func = obj_function.XFELTarget(mi=self.mi, dp=self.dp)
        self.objective_func.devices = []
        if not self.ui.cb_use_predef.checkState():
            self.objective_func.get_value = get_value_exp


    def set_m_status(self):

        state = self.is_le_addr_ok(self.ui.le_alarm)


        alarm_dev = str(self.ui.le_alarm.text())
        a_dev = AlarmDevice(alarm_dev)
        a_dev.mi = self.mi

        def is_ok():
            #alarm_dev = str(self.ui.le_alarm.text())
            alarm_min = self.ui.sb_alarm_min.value()
            alarm_max = self.ui.sb_alarm_max.value()
            #alarm_value = self.mi.get_value(alarm_dev)

            alarm_value = a_dev.get_value()

            print("ALARM: ", alarm_value, alarm_min, alarm_max)
            if alarm_min <= alarm_value <= alarm_max:
                return True
            return False
        self.m_status.is_ok = is_ok


    def getPlotData(self):
        """
        Collects data and updates plot on every GUI clock cycle.
        """
        #get times, penalties obj func data from the machine interface
        y = self.objective_func.penalties

        x = np.array(self.objective_func.times) - self.scanStartTime

        self.obj_func_line.setData(x=x, y=y)

        #plot data for all devices being scanned
        for dev in self.devices:
            y = np.array(dev.values)-self.multiPlotStarts[dev.eid]
            x = np.array(dev.time) - self.scanStartTime
            line = self.multilines[dev.eid]
            line.setData(x=x, y=y)


    def addPlots(self):
        """
        Initializes the GUIs plots and labels on startup.
        """
        #self.objective_func_pv = "test_obj"
        #setup plot 1 for obj func monitor
        self.plot1 = pg.PlotWidget(title="Objective Function Monitor", labels={'left': str(self.objective_func_pv), 'bottom':"Time (seconds)"})
        self.plot1.showGrid(1, 1, 1)
        self.plot1.getAxis('left').enableAutoSIPrefix(enable=False) # stop the auto unit scaling on y axes
        layout = QtGui.QGridLayout()
        self.ui.widget_2.setLayout(layout)
        layout.addWidget(self.plot1, 0, 0)

        #setup plot 2 for device monitor
        self.plot2 = pg.PlotWidget(title = "Device Monitor",labels={'left': "Device (Current - Start)", 'bottom': "Time (seconds)"})
        self.plot2.showGrid(1, 1, 1)
        self.plot2.getAxis('left').enableAutoSIPrefix(enable=False) # stop the auto unit scaling on y axes
        layout = QtGui.QGridLayout()
        self.ui.widget_3.setLayout(layout)
        layout.addWidget(self.plot2, 0, 0)

        #legend for plot 2
        self.leg2 = customLegend(offset=(75, 20))
        self.leg2.setParentItem(self.plot2.graphicsItem())

        #create the obj func line object
        color = QtGui.QColor(0, 255, 255)
        pen=pg.mkPen(color, width=3)
        self.obj_func_line = pg.PlotCurveItem(x=[], y=[], pen=pen, antialias=True)
        self.plot1.addItem(self.obj_func_line)

    def setUpMultiPlot(self, devices):
        """
        Reset plots when a new scan is started.
        """
        self.plot2.clear()
        self.multilines      = {}
        self.multiPvData     = {}
        self.multiPlotStarts = {}
        x = []
        y = []
        self.leg2.scene().removeItem(self.leg2)
        self.leg2 = customLegend(offset=(50, 10))
        self.leg2.setParentItem(self.plot2.graphicsItem())

        default_colors = [QtGui.QColor(255, 51, 51), QtGui.QColor(51, 255, 51), QtGui.QColor(255, 255, 51),QtGui.QColor(178, 102, 255)]
        for i, dev in enumerate(devices):

            #set the first 4 devices to have the same default colors
            if i < 4:
                color = default_colors[i]
            else:
                color = self.randColor()

            pen=pg.mkPen(color, width=2)
            self.multilines[dev.eid]  = pg.PlotCurveItem(x, y, pen=pen, antialias=True, name=str(dev.eid))
            self.multiPvData[dev.eid] = []
            self.multiPlotStarts[dev.eid] = dev.get_value()
            self.plot2.addItem(self.multilines[dev.eid])
            self.leg2.addItem(self.multilines[dev.eid], dev.eid, color=str(color.name()))

    def randColor(self):
        """
        Generate random line color for each device plotted.

        Returns:
                QColor object of a random color
        """
        hi = 255
        lo = 128
        c1 = np.random.randint(lo,hi)
        c2 = np.random.randint(lo,hi)
        c3 = np.random.randint(lo,hi)
        return QtGui.QColor(c1,c2,c3)

    def run_editor(self):
        #print(sys.modules[__name__].__file__, path)
        if platform.system() == 'Darwin':
            subprocess.call(['open', '-a', 'TextEdit', self.path_to_obj_func])
        elif platform.system() == 'Windows':
            subprocess.call(['C:\\Windows\\System32\\notepad.exe', self.path_to_obj_func])
        elif platform.system() == 'Linux':
            subprocess.call(['gedit', self.path_to_obj_func])
        else:
            print("Unknown platform")
            return
        self.set_obj_fun()

    def error_box(self, message):
        QtGui.QMessageBox.about(self, "Error box", message)



class customLegend(pg.LegendItem):
    """
    STUFF FOR PG CUSTOM LEGEND (subclassed from pyqtgraph).
    Class responsible for drawing a single item in a LegendItem (sans label).
    This may be subclassed to draw custom graphics in a Legend.
    """
    def __init__(self,size=None,offset=None):
        pg.LegendItem.__init__(self, size, offset)

    def addItem(self, item, name, color="CCFF00"):

        label = pg.LabelItem(name, color=color, size="6pt", bold=True)
        sample = None
        row = self.layout.rowCount()
        self.items.append((sample, label))
        self.layout.addItem(sample, row, 0)
        self.layout.addItem(label, row, 1)
        self.layout.setSpacing(0)


def main():

    """
    Funciton to start up the main program.

    Slecting a PV parameter set:
    If launched from the command line will take an argument with the filename of a parameter file.
    If no argv[1] is provided, the default list in ./parameters/lclsparams is used.

    Development mode:
    If devmode == False - GUI defaults to normal parameter list, defaults to nelder mead simplex
    if devmode == True  - GUI uses 4 development matlab PVs and loaded settings in the method "devmode()"
    """

    #try to get a pv list file name from commandline arg
    #this goes into initializing the reset panel PVs that show up in the GUI
    try:
        pvs = sys.argv[1]   # arg filename of params
    except:
        pvs = 'parameters/lcls_short.txt'#default filename

    #make pyqt threadsafe
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)

    #create the application
    app = QApplication(sys.argv)
    # setting the path variable for icon
    path = os.path.join(os.path.dirname(sys.modules[__name__].__file__), 'ocelot.png')
    app.setWindowIcon(QtGui.QIcon(path))
    window = OcelotInterfaceWindow()
    # setting the path variable for icon
    #path = os.path.join(os.path.dirname(sys.modules[__name__].__file__), 'ocelot.png')
    #window.setWindowIcon(QtGui.QIcon(path))

    #Build the PV list from dev PVs or selected source
    #window.ui.widget.getPvList(pvs)
    # set checkbot status
    #window.ui.widget.uncheckBoxes()

    timer = pg.QtCore.QTimer()
    timer.timeout.connect(window.scan_finished)
    timer.start(300)

    #show app

    #window.setWindowIcon(QtGui.QIcon('ocelot.png'))
    window.show()

    #Build documentaiton if source files have changed
    # TODO: make more universal
    #os.system("cd ./docs && xterm -T 'Ocelot Doc Builder' -e 'bash checkDocBuild.sh' &")
    #exit script
    sys.exit(app.exec_())

if __name__ == "__main__":

    main()