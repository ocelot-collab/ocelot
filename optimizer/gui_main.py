
from ocelot.optimizer.UIOcelotInterface_gen import *
import json
import scipy
from PyQt4.QtGui import QPixmap, QImage
from PIL import Image
import subprocess
import base64
from datetime import datetime

def send_to_desy_elog(author, title, severity, text, elog, image=None):
    """
    Send information to a supplied electronic logbook.
    Author Christopher Behrens (DESY)
    """

    # The DOOCS elog expects an XML string in a particular format. This string
    # is beeing generated in the following as an initial list of strings.
    succeded = True  # indicator for a completely successful job
    # list beginning
    elogXMLStringList = ['<?xml version="1.0" encoding="ISO-8859-1"?>',
                         '<entry>']
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
        #load in the dark theme style sheet
        self.loadStyleSheet()


        self.pb_save_as.clicked.connect(self.save_state_as)
        self.pb_load.clicked.connect(self.load_state_from)
        self.pb_rewrite.clicked.connect(self.rewrite_default)
        self.cb_use_isim.stateChanged.connect(self.change_state_scipy_setup)
        self.pb_hyper_file.clicked.connect(self.get_hyper_file)
        if scipy.__version__ < "0.18":
            self.g_box_isim.setEnabled(False)
            self.g_box_isim.setTitle("Initial Simplex does not work: scipy version: " + scipy.__version__)
            self.g_box_isim.setStyleSheet('QGroupBox  {color: red;}')

        self.le_a.textChanged.connect(self.check_address)

        self.sb_tdelay.valueChanged.connect(self.set_cycle)
        self.sb_ddelay.valueChanged.connect(self.set_cycle)
        #self.horizontalLayout_2.setStyleSheet("color: red")

            #font = self.pb_hyper_file.font()
        #font.setPointSize(16)
        #self.pb_hyper_file.setFont(font)
        #self.pb_hyper_file.setText("test")
        #self.pb_hyper_file.setStyleSheet("font: 16px, color: red")

        #self.window = window

    def set_cycle(self):
        """
        Select time for objective method data collection time.

        Scanner will wait this long to collect new data.
        """
        self.trim_delay = self.sb_tdelay.value()
        self.label_7.setText("Cycle Period = "+str(self.trim_delay ))
        self.Form.total_delay = self.trim_delay

    def check_address(self):
        self.is_le_addr_ok(self.le_a)
        self.is_le_addr_ok(self.le_b)
        self.is_le_addr_ok(self.le_c)
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
        #pvs = self.ui.widget.pvs
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
        table["max_iter"] = max_iter
        table["fun_a"] = fun_a
        table["fun_b"] = fun_b
        table["fun_c"] = fun_c
        table["obj_fun"] = obj_fun

        table["alarm_dev"] = alarm_dev
        table["alarm_min"] = alarm_min
        table["alarm_max"] = alarm_max

        table["seed_iter"] = self.sb_seed_iter.value()
        table["use_live_seed"] = self.cb_use_live_seed.checkState()

        table["isim_rel_step"] = self.sb_isim_rel_step.value()
        table["use_isim"] = self.cb_use_isim.checkState()

        table["hyper_file"] = self.Form.hyper_file

        table["set_best_sol"] = self.cb_set_best_sol.checkState()

        with open(filename, 'w') as f:
            json.dump(table, f)
        #pickle.dump(table, filename)
        print("SAVE State", table)

    def restore_state(self, filename):
        #try:
        with open(filename, 'r') as f:
            #data_new = pickle.load(f)
            table = json.load(f)

        # Build the PV list from dev PVs or selected source
        pvs = table["id"]
        print("pvs = ", pvs)
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
            self.sb_num_iter.setValue(max_iter)
            self.le_a.setText(fun_a)
            self.le_b.setText(fun_b)
            self.le_c.setText(fun_c)
            self.le_obf.setText(obj_fun)

            self.le_alarm.setText(table["alarm_dev"])
            self.sb_alarm_min.setValue(table["alarm_min"])
            self.sb_alarm_max.setValue(table["alarm_max"])

            self.sb_seed_iter.setValue(table["seed_iter"])
            self.cb_use_live_seed.setCheckState(table["use_live_seed"])

            self.sb_isim_rel_step.setValue(table["isim_rel_step"])
            self.cb_use_isim.setCheckState(table["use_isim"])
            self.change_state_scipy_setup()

            self.Form.hyper_file = table["hyper_file"]
            self.pb_hyper_file.setText(self.Form.hyper_file)

            self.cb_set_best_sol.setCheckState(table["set_best_sol"])
        except:
            pass


        print("RESTORE STATE")

    def save_state_as(self):
        filename = QtGui.QFileDialog.getSaveFileName(self.Form, 'Save State', filter ="txt (*.json *.)")
        if filename:
            self.Form.set_file = filename
            self.save_state(filename)

    def load_state_from(self):
        filename = QtGui.QFileDialog.getOpenFileName(self.Form, 'Load State', filter ="txt (*.json *.)")
        if filename:
            self.Form.set_file = filename
            self.restore_state(filename)

    def get_hyper_file(self):
        filename = QtGui.QFileDialog.getOpenFileName(self.Form, 'Load Hyper Parameters', filter ="txt (*.npy *.)")
        if filename:
            self.Form.hyper_file = str(filename)
            self.pb_hyper_file.setText(self.Form.hyper_file)
            #print(filename)


    def rewrite_default(self):
        self.Form.set_file = "default.json"
        self.save_state(self.Form.set_file)

    def logbook(self):
        filename = "screenshot"
        filetype = "png"
        self.screenShot(filename, filetype)
        table = self.Form.scan_params

        #curr_time = datetime.now()
        #timeString = curr_time.strftime("%Y-%m-%dT%H:%M:%S")
        text = ""
        if not self.cb_use_predef.checkState():
            text += "obj func: A   : " + str(self.le_a.text()).split("/")[-1] + "\n"
            text += "obj func: B   : " + str(self.le_b.text()).split("/")[-1] + "\n"
            text += "obj func: C   : " + str(self.le_c.text()).split("/")[-1] + "\n"
            text += "obj func: expr: " + str(self.le_obf.text()) + "\n"
        else:
            text += "obj func: A   : predefined  " + self.Form.objective_func.eid + "\n"
        if table != None:
            for i, dev in enumerate(table["devs"]):
                #print(dev.split("/"))
                text += "dev           : " + dev.split("/")[-1] + "   "+str(table["currents"][i][0]) +" --> " + str(table["currents"][i][1]) + "\n"

            text += "iterations    : " + str(table["iter"]) + "\n"
            text += "delay         : " + str(self.Form.total_delay) + "\n"
            text += "START-->STOP  :" + str(table["sase"][0]) +" --> " + str(table["sase"][1]) + "\n"
        #print(text)
        screenshot = open(self.Form.optimizer_path+filename+filetype)
        send_to_desy_elog(author="", title="OCELOT Optimization", severity="", text=text, elog="testlog", image=screenshot)


    def screenShot(self,filename,filetype):
        """
        Takes a screenshot of the whole gui window, saves png and ps images to file

        Args:
                fileName (str): Directory string of where to save the file
                filetype (str): String of the filetype to save
        """
        s = str(filename)+"."+str(filetype)
        p = QPixmap.grabWindow(self.Form.winId())
        p.save(s, 'png')
        #im = Image.open(s)
        #im.save(s[:-4]+".ps")
        p = p.scaled(465, 400)
        #save again a small image to use for the logbook thumbnail
        p.save(str(s[:-4])+"_sm.png", 'png')


    def loadStyleSheet(self):
        """ Sets the dark GUI theme from a css file."""
        try:
            self.cssfile = "style.css"
            with open(self.cssfile, "r") as f:
                self.Form.setStyleSheet(f.read())
        except IOError:
            print ('No style sheet found!')

    def change_state_scipy_setup(self):
        if self.cb_use_isim.checkState():
            self.label_23.setEnabled(True)
            self.sb_isim_rel_step.setEnabled(True)
        else:
            self.label_23.setEnabled(False)
            self.sb_isim_rel_step.setEnabled(False)

    def use_predef_fun(self):
        if self.cb_use_predef.checkState():
            self.le_a.setEnabled(False)
            self.le_b.setEnabled(False)
            self.le_c.setEnabled(False)
            self.le_obf.setEnabled(False)

            self.label_16.setEnabled(False)
            self.label_19.setEnabled(False)
            self.label_20.setEnabled(False)
            self.label_21.setEnabled(False)
        else:
            self.le_a.setEnabled(True)
            self.le_b.setEnabled(True)
            self.le_c.setEnabled(True)
            self.le_obf.setEnabled(True)

            self.label_16.setEnabled(True)
            self.label_19.setEnabled(True)
            self.label_20.setEnabled(True)
            self.label_21.setEnabled(True)

