"""
Most of GUI logic is placed here.
S.Tomin, 2017
"""

from ocelot.optimizer.UIOcelotInterface_gen import *
import json
import scipy
from PyQt5.QtGui import QPixmap, QImage
from PIL import Image
import subprocess
import base64
from datetime import datetime
import numpy as np
import sys
import os
import webbrowser
from shutil import copy


def send_to_desy_elog(author, title, severity, text, elog, image=None):
    """
    Send information to a supplied electronic logbook.
    Author: Christopher Behrens (DESY)
    """

    # The DOOCS elog expects an XML string in a particular format. This string
    # is beeing generated in the following as an initial list of strings.
    succeded = True  # indicator for a completely successful job
    # list beginning
    elogXMLStringList = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<entry>']

    # author information
    elogXMLStringList.append('<author>')
    elogXMLStringList.append(author)
    elogXMLStringList.append('</author>')
    # title information
    elogXMLStringList.append('<title>')
    elogXMLStringList.append(title)
    elogXMLStringList.append('</title>')
    # severity information
    elogXMLStringList.append('<severity>')
    elogXMLStringList.append(severity)
    elogXMLStringList.append('</severity>')
    # text information
    elogXMLStringList.append('<text>')
    elogXMLStringList.append(text)
    elogXMLStringList.append('</text>')
    # image information
    if image:
        try:
            encodedImage = base64.b64encode(image)
            elogXMLStringList.append('<image>')
            elogXMLStringList.append(encodedImage.decode())
            elogXMLStringList.append('</image>')
        except:  # make elog entry anyway, but return error (succeded = False)
            succeded = False
    # list end
    elogXMLStringList.append('</entry>')
    # join list to the final string
    elogXMLString = '\n'.join(elogXMLStringList)
    # open printer process
    try:
        lpr = subprocess.Popen(['/usr/bin/lp', '-o', 'raw', '-d', elog],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        # send printer job
        lpr.communicate(elogXMLString.encode('utf-8'))
    except:
        succeded = False
    return succeded


class MainWindow(Ui_Form):
    def __init__(self, Form):
        Ui_Form.__init__(self)
        self.setupUi(Form)
        self.Form = Form
        # load in the dark theme style sheet
        self.loadStyleSheet()
        self.widget.set_parent(Form)
        self.pb_save_as.clicked.connect(self.save_state_as)
        self.pb_load.clicked.connect(self.load_state_from)
        self.pb_rewrite.clicked.connect(self.rewrite_default)
        self.cb_use_isim.stateChanged.connect(self.change_state_scipy_setup)
        self.pb_hyper_file.clicked.connect(self.get_hyper_file)
        self.pb_logbook.clicked.connect(self.logbook)

        self.le_a.textChanged.connect(self.check_address)
        self.le_b.textChanged.connect(self.check_address)
        self.le_c.textChanged.connect(self.check_address)
        self.le_d.textChanged.connect(self.check_address)
        self.le_e.textChanged.connect(self.check_address)
        self.le_alarm.textChanged.connect(self.check_address)

        self.sb_tdelay.valueChanged.connect(self.set_cycle)
        self.sb_ddelay.valueChanged.connect(self.set_cycle)
        self.sb_nreadings.valueChanged.connect(self.set_cycle)
        self.cb_select_alg.currentIndexChanged.connect(self.change_state_scipy_setup)
        self.read_alarm = QtCore.QTimer()
        self.read_alarm.timeout.connect(self.alarm_value)
        self.read_alarm.start(1000)
        # self.horizontalLayout_2.setStyleSheet("color: red")

        # font = self.pb_hyper_file.font()
        # font.setPointSize(16)
        # self.pb_hyper_file.setFont(font)
        # self.pb_hyper_file.setText("test")
        # self.pb_hyper_file.setStyleSheet("font: 16px, color: red")

        # self.window = window

    def alarm_value(self):
        """
        reading alarm value
        :return:
        """
        dev = str(self.le_alarm.text())
        try:
            value = self.Form.mi.get_value(dev)
            self.label_alarm.setText(str(np.round(value, 2)))
        except:
            self.label_alarm.setText(str("None"))

    def set_cycle(self):
        """
        Select time for objective method data collection time.
        Scanner will wait this long to collect new data.
        """
        self.trim_delay = self.sb_tdelay.value()
        data_delay = self.sb_ddelay.value()*self.sb_nreadings.value()
        self.label_7.setText("Cycle Period = " + str(np.around(self.trim_delay + data_delay, 3)))
        self.Form.total_delay = self.trim_delay

    def check_address(self):
        self.is_le_addr_ok(self.le_a)
        self.is_le_addr_ok(self.le_b)
        self.is_le_addr_ok(self.le_c)
        self.is_le_addr_ok(self.le_d)
        self.is_le_addr_ok(self.le_e)
        self.is_le_addr_ok(self.le_alarm)

    def is_le_addr_ok(self, line_edit):
        dev = str(line_edit.text())
        state = True
        try:
            self.Form.mi.get_value(dev)
        except:
            state = False
        if state:
            line_edit.setStyleSheet("color: rgb(85, 255, 0);")
        else:
            line_edit.setStyleSheet("color: red")
        return state

    def save_state(self, filename):
        # pvs = self.ui.widget.pvs
        table = self.widget.get_state()

        table["use_predef"] = self.cb_use_predef.checkState()

        max_pen = self.sb_max_pen.value()
        timeout = self.sb_tdelay.value()


        max_iter = self.sb_num_iter.value()
        # objective function
        fun_a = str(self.le_a.text())
        fun_b = str(self.le_b.text())
        fun_c = str(self.le_c.text())
        obj_fun = str(self.le_obf.text())
        # alarm
        alarm_dev = str(self.le_alarm.text())
        alarm_min = self.sb_alarm_min.value()
        alarm_max = self.sb_alarm_max.value()

        table["max_pen"] = max_pen
        table["timeout"] = timeout
        table["nreadings"] = self.sb_nreadings.value()
        table["interval"] = self.sb_ddelay.value()

        table["max_iter"] = max_iter
        table["fun_a"] = fun_a
        table["fun_b"] = fun_b
        table["fun_c"] = fun_c
        table["fun_d"] = str(self.le_d.text())
        table["fun_e"] = str(self.le_e.text())
        table["obj_fun"] = obj_fun

        table["alarm_dev"] = alarm_dev
        table["alarm_min"] = alarm_min
        table["alarm_max"] = alarm_max
        table["alarm_timeout"] = self.sb_alarm_timeout.value()

        table["seed_iter"] = self.sb_seed_iter.value()
        table["use_live_seed"] = self.cb_use_live_seed.checkState()

        table["isim_rel_step"] = self.sb_isim_rel_step.value()
        table["use_isim"] = self.cb_use_isim.checkState()

        table["hyper_file"] = self.Form.hyper_file

        table["set_best_sol"] = self.cb_set_best_sol.checkState()

        table["algorithm"] = str(self.cb_select_alg.currentText())

        with open(filename, 'w') as f:
            json.dump(table, f)
        # pickle.dump(table, filename)
        print("SAVE State")

    def restore_state(self, filename):
        # try:
        with open(filename, 'r') as f:
            # data_new = pickle.load(f)
            table = json.load(f)

        # Build the PV list from dev PVs or selected source
        pvs = table["id"]
        self.widget.set_machine_interface(self.Form.mi, self.Form.dp)
        self.widget.getPvList(pvs)
        # set checkbot status
        self.widget.uncheckBoxes()
        self.widget.set_state(table)

        try:

            max_pen = table["max_pen"]
            timeout = table["timeout"]
            max_iter = table["max_iter"]
            fun_a = table["fun_a"]
            fun_b = table["fun_b"]
            fun_c = table["fun_c"]
            obj_fun = table["obj_fun"]

            if "use_predef" in table.keys(): self.cb_use_predef.setCheckState(table["use_predef"])
            self.sb_max_pen.setValue(max_pen)
            self.sb_tdelay.setValue(timeout)
            self.sb_nreadings.setValue(table["nreadings"])
            self.sb_ddelay.setValue(table["interval"])

            self.sb_num_iter.setValue(max_iter)
            self.le_a.setText(fun_a)
            self.le_b.setText(fun_b)
            self.le_c.setText(fun_c)
            self.le_d.setText(table["fun_d"])
            self.le_e.setText(table["fun_e"])
            self.le_obf.setText(obj_fun)

            self.le_alarm.setText(table["alarm_dev"])
            self.sb_alarm_min.setValue(table["alarm_min"])
            self.sb_alarm_max.setValue(table["alarm_max"])
            self.sb_alarm_timeout.setValue(table["alarm_timeout"])

            self.sb_seed_iter.setValue(table["seed_iter"])
            self.cb_use_live_seed.setCheckState(table["use_live_seed"])

            self.sb_isim_rel_step.setValue(table["isim_rel_step"])
            self.cb_use_isim.setCheckState(table["use_isim"])
            self.change_state_scipy_setup()

            self.Form.hyper_file = table["hyper_file"]
            self.pb_hyper_file.setText(self.Form.hyper_file)

            self.cb_set_best_sol.setCheckState(table["set_best_sol"])

            if "algorithm" in table.keys():
                index = self.cb_select_alg.findText(table["algorithm"], QtCore.Qt.MatchFixedString)

                if index >= 0:
                    self.cb_select_alg.setCurrentIndex(index)
            print("RESTORE STATE: OK")
        except:
            print("RESTORE STATE: ERROR")


    def save_state_as(self):

        filename = QtGui.QFileDialog.getSaveFileName(self.Form, 'Save State',
        self.Form.config_dir, "txt (*.json)", None, QtGui.QFileDialog.DontUseNativeDialog)[0]
        if filename:
            name = filename.split("/")[-1]
            parts = name.split(".")
            #print(parts)
            body_name = parts[0]

            if len(parts)<2 or parts[1] !="json":
                part = filename.split(".")[0]
                filename = part + ".json"
            copy(self.Form.obj_func_path, self.Form.obj_save_path + body_name +".py")
            #self.Form.set_file = filename
            self.save_state(filename)


    def load_state_from(self):
        filename = QtGui.QFileDialog.getOpenFileName(self.Form, 'Load State',
        self.Form.config_dir, "txt (*.json)", None, QtGui.QFileDialog.DontUseNativeDialog)[0]
        if filename:
            #print(filename)
            (body_name, extension) = filename.split("/")[-1].split(".")
            #print(self.Form.obj_save_path + body_name + ".py", self.Form.obj_func_path )
            copy(self.Form.obj_save_path + body_name + ".py", self.Form.obj_func_path )
            #self.Form.set_file = filename
            self.restore_state(filename)


    def get_hyper_file(self):
        #filename = QtGui.QFileDialog.getOpenFileName(self.Form, 'Load Hyper Parameters', filter="txt (*.npy *.)")
        filename = QtGui.QFileDialog.getOpenFileName(self.Form, 'Load Hyper Parameters',
        self.Form.optimizer_path  + "parameters", "txt (*.npy)", QtGui.QFileDialog.DontUseNativeDialog)
        if filename:
            self.Form.hyper_file = str(filename)
            self.pb_hyper_file.setText(self.Form.hyper_file)
            # print(filename)

    def rewrite_default(self):
        #self.Form.set_file = "default.json"
        self.save_state(self.Form.set_file)

    def logbook(self):
        """
        Method to send Optimization parameters + screenshot to eLogboob
        :return:
        """

        filename = "screenshot"
        filetype = "png"
        self.screenShot(filename, filetype)
        table = self.Form.scan_params

        # curr_time = datetime.now()
        # timeString = curr_time.strftime("%Y-%m-%dT%H:%M:%S")
        text = ""

        if not self.cb_use_predef.checkState():
            text += "obj func: A   : " + str(self.le_a.text()).split("/")[-2]  + "/"+ str(self.le_a.text()).split("/")[-1] + "\n"
            if str(self.le_b.text()) != "":
                text += "obj func: B   : " + str(self.le_b.text()).split("/")[-2] + "/" + str(self.le_b.text()).split("/")[-1] + "\n"
            if str(self.le_c.text()) != "":
                text += "obj func: C   : " + str(self.le_c.text()).split("/")[-2] + "/" + str(self.le_c.text()).split("/")[-1] + "\n"
            if str(self.le_d.text()) != "":
                text += "obj func: D   : " + str(self.le_d.text()).split("/")[-2] + "/" + str(self.le_d.text()).split("/")[-1] + "\n"
            if str(self.le_e.text()) != "":
                text += "obj func: E   : " + str(self.le_e.text()).split("/")[-2] + "/" + str(self.le_e.text()).split("/")[-1] + "\n"
            text += "obj func: expr: " + str(self.le_obf.text()) + "\n"
        else:
            text += "obj func: A   : predefined  " + self.Form.objective_func.eid + "\n"
        if table != None:
            for i, dev in enumerate(table["devs"]):
                # print(dev.split("/"))
                text += "dev           : " + dev.split("/")[-2] + "/" + dev.split("/")[-1] + "   " + str(table["currents"][i][0]) + " --> " + str(
                    table["currents"][i][1]) + "\n"

            text += "iterations    : " + str(table["iter"]) + "\n"
            text += "delay         : " + str(self.Form.total_delay) + "\n"
            text += "START-->STOP  : " + str(table["sase"][0]) + " --> " + str(table["sase"][1]) + "\n"
            text += "Method        : " + str(table["method"]) + "\n"
        #print("table", table)
        #print(text)
        screenshot = open(self.Form.optimizer_path + filename + "." + filetype, 'rb')
        #print(screenshot)
        res = send_to_desy_elog(author="", title="OCELOT Optimization", severity="INFO", text=text, elog=self.Form.logbook,
                          image=screenshot.read())

        if not res:
            self.Form.error_box("error during eLogBook sending")

    def screenShot(self, filename, filetype):

        """
        Takes a screenshot of the whole gui window, saves png and ps images to file
        :param filename: (str) Directory string of where to save the file
        :param filetype: (str) String of the filetype to save
        :return:
        """

        s = str(filename) + "." + str(filetype)
        p = QPixmap.grabWindow(self.Form.winId())
        p.save(s, 'png')
        # im = Image.open(s)
        # im.save(s[:-4]+".ps")
        p = p.scaled(465, 400)
        # save again a small image to use for the logbook thumbnail
        p.save(str(s[:-4]) + "_sm.png", 'png')

    def loadStyleSheet(self):
        """
        Sets the dark GUI theme from a css file.
        :return:
        """
        try:
            self.cssfile = "style.css"
            with open(self.cssfile, "r") as f:
                self.Form.setStyleSheet(f.read())
        except IOError:
            print ('No style sheet found!')

    def change_state_scipy_setup(self):
        """
        Method to enable/disable "Scipy Scanner Setup". If scipy version < "0.18" then QGroup will be disable.
        :return:
        """
        #print("SCIPY", str(self.cb_select_alg.currentText()))
        if scipy.__version__ < "0.18" and str(self.cb_select_alg.currentText()) == self.Form.name_simplex:
            #self.cb_use_isim.setCheckState(False)
            self.g_box_isim.setEnabled(False)
            self.g_box_isim.setTitle("Initial Simplex does not work: scipy version: " + scipy.__version__)
            self.g_box_isim.setStyleSheet('QGroupBox  {color: red;}')
        elif scipy.__version__ >= "0.18" and str(self.cb_select_alg.currentText()) == self.Form.name_simplex:
            #print(str(self.cb_select_alg.currentText()))
            self.g_box_isim.setEnabled(True)
            self.cb_use_isim.setEnabled(True)
            self.g_box_isim.setTitle("Simplex/Scipy Scanner Setup")
            self.g_box_isim.setStyleSheet('QGroupBox  {color: white;}')

        if self.cb_use_isim.checkState():
            self.label_23.setEnabled(True)
            self.sb_isim_rel_step.setEnabled(True)
        else:
            self.label_23.setEnabled(False)
            self.sb_isim_rel_step.setEnabled(False)

        if str(self.cb_select_alg.currentText()) == self.Form.name_custom:
            self.g_box_isim.setEnabled(True)
            self.label_23.setEnabled(True)
            self.sb_isim_rel_step.setEnabled(True)
            self.g_box_isim.setTitle("Custom Minimizer Scanner Setup")
            self.g_box_isim.setStyleSheet('QGroupBox  {color: white;}')
            #self.cb_use_isim.setCheckState(True)
            self.cb_use_isim.setEnabled(False)
            self.sb_isim_rel_step.setValue(5)

        if str(self.cb_select_alg.currentText()) in [self.Form.name_simplex_norm, self.Form.name_gauss_sklearn]:
            self.g_box_isim.setEnabled(True)
            self.label_23.setEnabled(True)
            self.sb_isim_rel_step.setEnabled(True)
            self.g_box_isim.setTitle("Simplex With Normalization")
            self.g_box_isim.setStyleSheet('QGroupBox  {color: white;}')
            #self.cb_use_isim.setCheckState(True)
            self.cb_use_isim.setEnabled(False)
            self.sb_isim_rel_step.setValue(5)


    def use_predef_fun(self):
        if self.cb_use_predef.checkState():
            self.le_a.setEnabled(False)
            self.le_b.setEnabled(False)
            self.le_c.setEnabled(False)
            self.le_d.setEnabled(False)
            self.le_e.setEnabled(False)
            self.le_obf.setEnabled(False)

            self.label_16.setEnabled(False)
            self.label_19.setEnabled(False)
            self.label_20.setEnabled(False)
            self.label_21.setEnabled(False)
            self.label_28.setEnabled(False)
            self.label_29.setEnabled(False)
        else:
            self.le_a.setEnabled(True)
            self.le_b.setEnabled(True)
            self.le_c.setEnabled(True)
            self.le_d.setEnabled(True)
            self.le_e.setEnabled(True)
            self.le_obf.setEnabled(True)

            self.label_16.setEnabled(True)
            self.label_19.setEnabled(True)
            self.label_20.setEnabled(True)
            self.label_21.setEnabled(True)
            self.label_28.setEnabled(True)
            self.label_29.setEnabled(True)

    def open_help(self):
        """
        method to open the Help in the webbrowser
        :return: None
        """

        if sys.platform == 'win32':
            url = self.Form.optimizer_path+"docs\\_build\\html\\index.html"
            #os.startfile(url)
            webbrowser.open(url)
        elif sys.platform == 'darwin':
            url = "file://"+self.Form.optimizer_path+"docs/_build/html/index.html"
            webbrowser.open(url)
            #subprocess.Popen(['open', url])
        else:
            url = "file://" + self.Form.optimizer_path + "docs/_build/html/index.html"
            try:
                subprocess.Popen(['xdg-open', url])
            except OSError:
                print('Please open a browser on: ' + url)

