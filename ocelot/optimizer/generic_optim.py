#!/opt/anaconda4/bin/python
"""
This is deep modification of SLAC version of the Ocelot GUI for the European XFEL facility.
Sergey Tomin, 2017.
Ocelot GUI, interface for running and testing accelerator optimization methods
This file primarily contains the code for the UI and GUI
The scanner classes are contained in an external file, scannerThreads.py
The resetpanel widget is also contained in a separate module, resetpanel
Tyler Cope, 2016
"""
from __future__ import absolute_import, print_function
import json
import sys
import os
import sklearn
sklearn_version = sklearn.__version__

path = os.path.realpath(__file__)
indx = path.find("ocelot/optimizer")
print("PATH", os.path.realpath(__file__))
sys.path.append(path[:indx])

# for pyqtgraph import
#sys.path.append(path[:indx]+"ocelot")



from PyQt5.QtGui import QPixmap
from PyQt5 import QtGui, QtCore
from PyQt5.QtWidgets import QApplication, QFrame, QMessageBox, QMainWindow, QDialog
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


from ocelot.optimizer.gui_main import *

from ocelot.optimizer.mint.opt_objects import *
from ocelot.optimizer.mint import mint
from ocelot.optimizer.mint import opt_objects as obj

from ocelot.utils import db
from ocelot.optimizer.mint.xfel_interface import *


class OcelotInterfaceWindow(QFrame):
    """ Main class for the GUI application """
    def __init__(self):
        """
        Initialize the GUI and QT UI aspects of the application.
        Initialize the scan parameters.
        Connect start and logbook buttons on the scan panel.
        Initialize the plotting.
        Make the timer object that updates GUI on clock cycle during a scan.
        """
        # PATHS
        path = os.path.realpath(__file__)
        indx = path.find("ocelot" + os.sep + "optimizer")
        self.path2ocelot = path[:indx]
        self.optimizer_path = self.path2ocelot + "ocelot" + os.sep + "optimizer" + os.sep
        self.config_dir = self.path2ocelot + "config_optim" +os.sep
        self.set_file = self.config_dir + "default.json" # ./parameters/default.json"
        self.obj_func_path = self.optimizer_path + "mint" + os.sep + "obj_function.py"
        self.obj_save_path = self.config_dir +  "obj_funcs" + os.sep
        # initialize
        QFrame.__init__(self)

        self.logbook = "xfellog"
        self.dev_mode = False

        self.ui = MainWindow(self)

        self.name_simplex = "Nelder-Mead Simplex"
        self.name_gauss = "Gaussian Process"
        self.name_gauss_sklearn = "Gaussian Process sklearn"
        self.name_custom = "Custom Minimizer"
        self.name_simplex_norm = "Simplex Norm."
        # self.name4 = "Conjugate Gradient"
        # self.name5 = "Powell's Method"
        # switch of GP and custom Mininimizer
        self.ui.cb_select_alg.addItem(self.name_simplex)
        #self.ui.cb_select_alg.addItem(self.name_gauss)
        self.ui.cb_select_alg.addItem(self.name_custom)
        self.ui.cb_select_alg.addItem(self.name_simplex_norm)
        if sklearn_version >= "0.18":
            self.ui.cb_select_alg.addItem(self.name_gauss_sklearn)


        #self.ui.pb_help.clicked.connect(lambda: os.system("firefox file://"+self.optimizer_path+"docs/build/html/index.html"))
        self.ui.pb_help.clicked.connect(self.ui.open_help)
        if self.dev_mode:
            self.mi = TestMachineInterface()
        else:
            self.mi = XFELMachineInterface()
        self.dp = TestDeviceProperties(ui=self.ui.widget)

        self.total_delay = self.ui.sb_tdelay.value()
        self.opt_control = mint.OptControl()
        self.opt_control.alarm_timeout = self.ui.sb_alarm_timeout.value()

        #self.objective_func = obj_function.XFELTarget()
        self.objective_func_pv = "test_obj"

        self.show_obj_value = False
        self.addPlots()

        # database
        self.dbname =  self.config_dir + "test.db" #"./parameters/test.db"
        # db.create_db(self.dbname)
        print(self.optimizer_path, self.set_file, self.dbname)
        try:
            self.db = db.PerfDB(dbname=self.dbname)
        except:
            self.db = None
            self.error_box(message="Database is not available")

        self.scan_params = None
        self.hyper_file = "../parameters/hyperparameters.npy"

        self.ui.pb_start_scan.clicked.connect(self.start_scan)
        self.ui.pb_edit_obj_func.clicked.connect(self.run_editor)
        self.ui.cb_use_predef.stateChanged.connect(self.set_obj_fun)

        self.ui.restore_state(self.set_file)
        self.path_to_obj_func = os.path.join(os.path.dirname(sys.modules[__name__].__file__), 'mint/obj_function.py')


        self.set_obj_fun()
        self.m_status = mint.MachineStatus()
        self.set_m_status()
        self.opt_control = mint.OptControl()
        self.opt_control.m_status = self.m_status


        #timer for plots, starts when scan starts
        self.multiPvTimer = QtCore.QTimer()
        self.multiPvTimer.timeout.connect(self.getPlotData)


    def scan_method_select(self):
        """
        Sets scanner method from options panel combo box selection.
        This method executes from the runScan() method, when the UI "Start Scan" button is pressed.
        :return: Selected scanner object
                 These objects are contrained in the scannerThreads.py file
        """
        current_method = self.ui.cb_select_alg.currentText()

        #GP Method
        if current_method == self.name_gauss:
            minimizer = mint.GaussProcess()

        elif current_method == self.name_gauss_sklearn:
            minimizer = mint.GaussProcessSKLearn()
        # Custom Minimizer
        elif current_method == self.name_custom:
            minimizer = mint.CustomMinimizer()

        elif current_method == self.name_simplex_norm:
            minimizer = mint.Simplex()

        #simplex Method
        else:
            minimizer = mint.Simplex()

        self.method_name = minimizer.__class__.__name__

        return minimizer


    def closeEvent(self, event):
        self.ui.save_state(self.set_file)
        if self.ui.pb_start_scan.text() == "Stop optimization":
            self.opt.opt_ctrl.stop()
            self.m_status.is_ok = lambda: True
            del(self.opt)
            self.ui.pb_start_scan.setStyleSheet("color: rgb(85, 255, 127);")
            self.ui.pb_start_scan.setText("Start optimization")
            return 0
        QFrame.closeEvent(self, event)

    def start_scan(self):
        """
        Method to start/stop the Optimizer.
        """
        self.scanStartTime = time.time()

        if self.ui.pb_start_scan.text() == "Stop optimization":
            # stop the optimization
            self.opt.opt_ctrl.stop()

            self.m_status.is_ok = lambda: True
            # Save the optimization parameters to the database
            self.save2db()
            del(self.opt)
            # Setting the button
            self.ui.pb_start_scan.setStyleSheet("color: rgb(85, 255, 127);")
            self.ui.pb_start_scan.setText("Start optimization")
            return 0

        self.pvs = self.ui.widget.getPvsFromCbState()
        self.devices = self.ui.widget.get_devices(self.pvs)

        if len(self.devices) == 0:
            self.error_box(message="Check Devices")
            return 0
        for dev in self.devices:
            val = dev.get_value()
            if dev.check_limits(val):
                self.error_box(message="Check the Limits")
                return 0
        self.setUpMultiPlot(self.devices)
        self.multiPvTimer.start(100)

        # set the Objective function from GUI or from file mint.obj_function.py (reloading)
        self.set_obj_fun()

        # Set minimizer - the optimization method (Simplex, GP, ...)
        minimizer = self.scan_method_select()

        # configure the Minimizer
        if minimizer.__class__ in [mint.GaussProcess, mint.GaussProcessSKLearn]:
            minimizer.seed_iter = self.ui.sb_seed_iter.value()
            minimizer.seed_timeout = self.ui.sb_tdelay.value()
            minimizer.hyper_file = self.hyper_file
            minimizer.norm_coef = self.ui.sb_isim_rel_step.value()/ 100.

        elif minimizer.__class__ == mint.Simplex:
            if self.ui.cb_use_isim.checkState():
                minimizer.dev_steps = []

                for dev in self.devices:
                    if dev.simplex_step == 0:
                        lims = dev.get_limits()
                        rel_step = self.ui.sb_isim_rel_step.value()
                        minimizer.dev_steps.append((lims[1] - lims[0]) * rel_step / 100.)
            else:
                minimizer.dev_steps = None

        elif minimizer.__class__ == mint.CustomMinimizer:
            minimizer.dev_steps = []

            for dev in self.devices:
                if dev.simplex_step == 0:
                    lims = dev.get_limits()
                    rel_step = self.ui.sb_isim_rel_step.value()
                    print(dev.id, rel_step)
                    minimizer.dev_steps.append((lims[1] - lims[0]) * rel_step / 100.)
            print("MINImizer steps", minimizer.dev_steps)

        self.max_iter = self.ui.sb_num_iter.value()
        minimizer.max_iter = self.max_iter

        # Optimizer initialization
        self.opt = mint.Optimizer()

        if self.ui.cb_select_alg.currentText() == self.name_simplex_norm:
            self.opt.normalization = True
            self.opt.norm_coef = self.ui.sb_isim_rel_step.value()*0.01
            print("OPT", self.opt.norm_coef)
        # Option - set best solution after optimization or not
        self.opt.set_best_solution = self.ui.cb_set_best_sol.checkState()

        self.set_m_status()
        self.opt_control.m_status = self.m_status
        self.opt.opt_ctrl = self.opt_control
        self.opt.timeout = self.total_delay

        self.opt.minimizer = minimizer

        seq = [mint.Action(func=self.opt.max_target_func, args=[self.objective_func, self.devices])]
        self.opt.seq = seq

        #self.opt.eval(seq)
        self.opt.start()

        # Setting the button
        self.ui.pb_start_scan.setText("Stop optimization")
        self.ui.pb_start_scan.setStyleSheet("color: red")

    def scan_finished(self):
        try:
            if self.ui.pb_start_scan.text() == "Stop optimization" and not (self.opt.isAlive()):
                self.ui.pb_start_scan.setStyleSheet("color: rgb(85, 255, 127);")
                self.ui.pb_start_scan.setText("Start optimization")
                self.save2db()
                print("scan_finished: OK")
        except:
            print("scan_finished: ERROR")

    def create_devices(self, pvs):
        """
        Method to create devices using only channels (PVs)

        :param pvs: str, device address/channel/PV
        :return: list of the devices [mint.opt_objects.Device(eid=pv[0]), mint.opt_objects.Device(eid=pv[1]), ... ]
        """
        # TODO: add new method for creation of devices
        devices = []
        for pv in pvs:
            if self.dev_mode:
                dev = obj.TestDevice(eid=pv)
            else:
                dev = obj.Device(eid=pv)
            dev.mi = self.mi
            dev.dp = self.dp
            devices.append(dev)
        return devices


    def indicate_machine_state(self):
        """
        Method to indicate of the machine status. Red frames around graphics means that machine status is not OK.
        :return:
        """
        if not self.opt_control.is_ok:
            if self.ui.widget_3.styleSheet() == "background-color:red;":
                return
            self.ui.widget_2.setStyleSheet("background-color:red;")
            self.ui.widget_3.setStyleSheet("background-color:red;")

        else:
            if self.ui.widget_3.styleSheet() == "background-color:323232;":
                return
            self.ui.widget_2.setStyleSheet("background-color:323232;")
            self.ui.widget_3.setStyleSheet("background-color:323232;")


    def save2db(self):
        """
        Save optimization parameters to the Database

        :return: None
        """
        dump2json = {}
        self.scan_params = {"devs": [], "currents": [], "iter":0, "sase": [0,0],"pen":[0,0], "obj":[]}
        d_names = []
        d_start = []
        d_stop = []
        for dev in self.devices:
            self.scan_params["devs"].append(dev.eid)
            d_names.append(dev.eid + "_val")
            d_start.append(dev.values[0])
            d_stop.append(dev.values[-1])

            self.scan_params["currents"].append([dev.values[0], dev.values[-1]])

            d_names.append(dev.eid + "_lim")
            d_start.append(dev.get_limits()[0])
            d_stop.append(dev.get_limits()[1])
            dump2json[dev.eid] = dev.values

        self.scan_params["iter"] = len(self.objective_func.penalties)

        self.scan_params["sase"] = [self.objective_func.values[0], self.objective_func.values[-1]]
        self.scan_params["pen"] = [self.objective_func.penalties[0], self.objective_func.penalties[-1]]
        self.scan_params["method"] = self.method_name

        o_names = ["obj_id", "obj_value", "obj_pen", "niter"]
        o_start = [self.objective_func.eid, self.objective_func.values[0], self.objective_func.penalties[0], self.max_iter]
        o_stop = [self.objective_func.eid, self.objective_func.values[-1], self.objective_func.penalties[-1],  len(self.objective_func.penalties)]


        start_sase = self.objective_func.values[0]
        stop_sase = self.objective_func.values[-1]

        #print('current actions in tuning', [(t.id, t.tuning_id, t.sase_start, t.sase_end) for t in self.db.get_actions()])


        # add new data here START
        new_data_name =  ["method"]
        new_data_start = [self.method_name]
        new_data_end =   [self.method_name]
        # add new data here END

        param_names = o_names + d_names + new_data_name
        start_vals = o_start+d_start + new_data_start
        end_vals = o_stop+d_stop + new_data_end

        try:
            self.db.new_tuning({'wl': 13.6, 'charge': self.mi.get_charge(), 'comment': 'test tuning'})
            tune_id = self.db.current_tuning_id()
            self.db.new_action(tune_id, start_sase=start_sase, end_sase=stop_sase)

            action_id = self.db.current_action_id()
            self.db.add_action_parameters(tune_id, action_id, param_names=param_names, start_vals=start_vals,
                                     end_vals=end_vals)

            print('current actions', [(t.id, t.tuning_id, t.sase_start, t.sase_end) for t in self.db.get_actions()])
            print('current action parameters', [(p.tuning_id, p.action_id, p.par_name, p.start_value, p.end_value) for p in
                                                self.db.get_action_parameters(tune_id, action_id)])
        except:
            self.error_box(message="Database error")

        dump2json["dev_times"] = self.devices[0].times
        dump2json["obj_values"] = self.objective_func.values
        dump2json["obj_times"] = self.objective_func.times

        #path = os.getcwd()
        #indx = path.find("ocelot")
        #path = path[:indx]
        path = self.path2ocelot
        print("JSON", path)
        filename = os.path.join(path, "data", datetime.now().strftime("%Y-%m-%d %H-%M-%S") + ".json")
        #print(filename)
        try:
            with open(filename, 'w+') as f:
                json.dump(dump2json, f)
        except:
            print("ERROR. Could not write history")

    def set_obj_fun(self):
        """
        Method to set objective function from the GUI (channels A,B,C) or reload module obj_function.py

        :return: None
        """
        try:
            from ocelot.optimizer.mint import obj_function
            reload(obj_function)
            self.ui.pb_edit_obj_func.setStyleSheet("background: #4e4e4e")
        except:
            self.ui.pb_edit_obj_func.setStyleSheet("background: red")
            self.ui.cb_use_predef.setCheckState(False)
            print("ERROR set objective function")

        self.ui.use_predef_fun()

        if self.ui.cb_use_predef.checkState():
            print("RELOAD Module Objective Function")
            self.objective_func = obj_function.XFELTarget(mi=self.mi, dp=self.dp)
            self.objective_func.devices = []
        else:
            # disable button "Edit Objective Function"
            # self.ui.pb_edit_obj_func.setEnabled(False)
            line_edits = [self.ui.le_a, self.ui.le_b, self.ui.le_c, self.ui.le_d, self.ui.le_e]

            a_str = str(self.ui.le_a.text())
            state_a = self.ui.is_le_addr_ok(self.ui.le_a)

            b_str = str(self.ui.le_b.text())
            state_b = self.ui.is_le_addr_ok(self.ui.le_b)

            c_str = str(self.ui.le_c.text())
            state_c = self.ui.is_le_addr_ok(self.ui.le_c)

            d_str = str(self.ui.le_d.text())
            state_d = self.ui.is_le_addr_ok(self.ui.le_d)

            e_str = str(self.ui.le_e.text())
            state_e = self.ui.is_le_addr_ok(self.ui.le_e)

            func = str(self.ui.le_obf.text())

            def get_value_exp():
                A = 0.
                B = 0.
                C = 0.
                D = 0.
                E = 0.
                if state_a:
                    A = self.mi.get_value(a_str)
                if state_b:
                    B = self.mi.get_value(b_str)
                if state_c:
                    C = self.mi.get_value(c_str)
                if state_d:
                    D = self.mi.get_value(d_str)
                if state_e:
                    E = self.mi.get_value(e_str)
                return eval(func)

            self.objective_func = obj.Target(eid=a_str)
            self.objective_func.devices = []
            self.objective_func.get_value = get_value_exp

        # set maximum penalty
        self.objective_func.pen_max = self.ui.sb_max_pen.value()
        # set number of the readings
        self.objective_func.nreadings = self.ui.sb_nreadings.value()
        # set interval between readings
        self.objective_func.interval = self.ui.sb_ddelay.value()
        if self.dev_mode:
            def get_value_dev_mode():
                values = np.array([dev.get_value() for dev in self.devices])
                print("I am here!")
                return np.sum(np.exp(-np.power((values - np.ones_like(values)), 2) / 5.))

            self.objective_func.get_value = get_value_dev_mode

    def set_m_status(self):
        """
        Method to set the MachineStatus method self.is_ok using GUI Alarm channel and limits
        :return: None
        """

        alarm_dev = str(self.ui.le_alarm.text()).replace(" ", "")
        print("alarm_dev", alarm_dev)
        if alarm_dev == "":
            return

        state = self.ui.is_le_addr_ok(self.ui.le_alarm)
        #state = False

        a_dev = AlarmDevice(alarm_dev)
        a_dev.mi = self.mi
        print(a_dev)
        if not state:
            def is_ok():
                print("ALARM switched off")
                return True
        else:
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
        if len(self.objective_func.times) == 0:
            return
        y = self.objective_func.penalties

        x = np.array(self.objective_func.times) - self.objective_func.times[0]

        self.obj_func_line.setData(x=x, y=y)
        if self.show_obj_value:
            self.obj_func_value.setData(x=x, y=self.objective_func.values)

        #plot data for all devices being scanned
        for dev in self.devices:
            if len(dev.times) == 0:
                return
            y = np.array(dev.values)-self.multiPlotStarts[dev.eid]
            x = np.array(dev.times) - np.array(dev.times)[0]
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
        self.plot2 = pg.PlotWidget(title="Device Monitor", labels={'left': "Device (Current - Start)", 'bottom': "Time (seconds)"})
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
        self.obj_func_value = pg.PlotCurveItem(x=[], y=[], pen=pg.mkPen(QtGui.QColor(255, 255, 51), width=3),
                                               antialias=True, name="value")
        self.plot1.addItem(self.obj_func_line)
        if self.show_obj_value:
            self.plot1.addItem(self.obj_func_value)

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
            self.multilines[dev.eid] = pg.PlotCurveItem(x, y, pen=pen, antialias=True, name=str(dev.eid))
            self.multiPvData[dev.eid] = []
            self.multiPlotStarts[dev.eid] = dev.get_value()
            self.plot2.addItem(self.multilines[dev.eid])
            self.leg2.addItem(self.multilines[dev.eid], dev.eid, color=str(color.name()))

    def randColor(self):
        """
        Generate random line color for each device plotted.
        :return: QColor object of a random color
        """
        hi = 255
        lo = 128
        c1 = np.random.randint(lo,hi)
        c2 = np.random.randint(lo,hi)
        c3 = np.random.randint(lo,hi)
        return QtGui.QColor(c1,c2,c3)

    def run_editor(self):
        """
        Run the editor for edition of the objective function in obj_function.py
        :return:
        """
        if platform.system() == 'Darwin':
            #subprocess.call(['open', '-a', 'TextEdit', self.path_to_obj_func])
            subprocess.call(['open', self.path_to_obj_func])
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
        #QtGui.QMessageBox.critical(self, "Error box", message)

class customLegend(pg.LegendItem):
    """
    STUFF FOR PG CUSTOM LEGEND (subclassed from pyqtgraph).
    Class responsible for drawing a single item in a LegendItem (sans label).
    This may be subclassed to draw custom graphics in a Legend.
    """
    def __init__(self, size=None, offset=None):
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
    #try:
    #    pvs = sys.argv[1]   # arg filename of params
    #except:
    #    pvs = 'parameters/lcls_short.txt'#default filename

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

    indicator = QtCore.QTimer()
    indicator.timeout.connect(window.indicate_machine_state)
    indicator.start(10)

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
