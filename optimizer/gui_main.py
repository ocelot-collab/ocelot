
from ocelot.optimizer.UIOcelotInterface_gen import *
import json
from PyQt4.QtGui import QPixmap, QImage
from PIL import Image

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
        #font = self.pb_hyper_file.font()
        #font.setPointSize(16)
        #self.pb_hyper_file.setFont(font)
        #self.pb_hyper_file.setText("test")
        #self.pb_hyper_file.setStyleSheet("font: 16px, color: red")

        #self.window = window

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

            self.Form.hyper_file = table["hyper_file"]
            self.pb_hyper_file.setText(self.Form.hyper_file)

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
        self.screenShot( filename, filetype)

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
        p = p.scaled(465,400)
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

