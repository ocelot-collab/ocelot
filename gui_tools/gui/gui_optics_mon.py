from gui.UIOpticsMonitor import *
#from table_test import *
import sys
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QAction, QTableWidget, QTableWidgetItem, QVBoxLayout
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from PyQt5 import QtCore
from gui.device.device_ui import *
from ocelot import *

class OcelotOpticsWindow(QMainWindow):
    """ Main class for the GUI application """
    def __init__(self, master=None):
        super().__init__()

        # initialize
        #QFrame.__init__(self)

        #self.logbook = "xfellog"
        #self.dev_mode = False
        self.master = master
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.ui.menuBar.setNativeMenuBar(False)
        self.ui.mainToolBar.setVisible(False)
        self.loadStyleSheet()
        self.ui.gridLayout_6.setContentsMargins(0,0,0,0)
        self.ui.pb_hide_show.clicked.connect(self.hide_show_quads)

        self.ui.pb_start_stop.clicked.connect(self.start_stop_reading)

        self.cyclic_read_timer = QtCore.QTimer()
        self.cyclic_read_timer.timeout.connect(self.master.auto_reading)

        self.ui.w_q_table.ui.pb_calculate.clicked.connect(self.master.calc_twiss)
        self.ui.cb_design.setChecked(True)
        self.ui.cb_design.stateChanged.connect(self.check_boxes_design)
        self.ui.cb_otrc55.stateChanged.connect(self.check_boxes)
        self.ui.cb_otrb218.stateChanged.connect(self.check_boxes)
        self.ui.cb_otrb450.stateChanged.connect(self.check_boxes)

    def check_boxes(self):
        if self.ui.cb_otrb218.isChecked():
            self.ui.cb_design.setChecked(False)
        elif self.ui.cb_otrb450.isChecked():
            self.ui.cb_design.setChecked(False)
        else:
            self.ui.cb_design.setChecked(False)

    def check_boxes_design(self):
        self.ui.cb_otrc55.setChecked(False)
        self.ui.cb_otrb218.setChecked(False)
        self.ui.cb_otrb450.setChecked(False)

    def stop_cyclic_reading(self):
        self.cyclic_read_timer.stop()
        self.ui.pb_start_stop.setStyleSheet("color: rgb(85, 255, 127);")
        self.ui.pb_start_stop.setText("Start Cyclic Reading")

    def start_stop_reading(self):
        """
        Method to start/stop feedback timer.
        sb_feedback_sec - spinBox - set seconds for timer
        pb_feedback - pushBatton Off/On
        feedback_timer - timer
        :return:
        """
        # print("I am here")

        # if self.ui.pb_start_statistics.text() == "Statistics Accum On":
        #    return 0

        # self.orbit = self.orbit_class.create_Orbit_obj()

        delay = self.ui.sb_delay.value() * 1000
        if self.ui.pb_start_stop.text() == "Stop Cyclic Reading":
            self.stop_cyclic_reading()
        else:
            self.cyclic_read_timer.start(delay)
            self.ui.pb_start_stop.setText("Stop Cyclic Reading")
            self.ui.pb_start_stop.setStyleSheet("color: red")


    def hide_show_quads(self):

        if self.ui.pb_hide_show.text() == "Hide Table":
            #self.ui.w_q_table.setParent(None)  # remove old table
            self.ui.w_q_table.hide()
            self.ui.pb_hide_show.setText("Show Table")
            self.ui.pb_hide_show.setStyleSheet("color: rgb(255, 0, 255);")
        else:
            self.ui.w_q_table.show()
            #self.ui.w_q_table = VOclTable(self.ui.tab_4)
            #self.ui.w_q_table.setEnabled(True)
            #sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
            #sizePolicy.setHorizontalStretch(0)
            #sizePolicy.setVerticalStretch(0)
            #sizePolicy.setHeightForWidth(self.ui.w_q_table.sizePolicy().hasHeightForWidth())
            #self.ui.w_q_table.setSizePolicy(sizePolicy)
            #self.ui.w_q_table.setMinimumSize(QtCore.QSize(420, 200))
            #self.ui.w_q_table.setMaximumSize(QtCore.QSize(420, 16777215))
            #self.ui.gridLayout_6.addWidget(self.ui.w_q_table, 0, 3, 1, 1)
            #self.init_quad_table(self.master.quads, calc_obj=self.master.calc_twiss)
            self.ui.pb_hide_show.setText("Hide Table")
            self.ui.pb_hide_show.setStyleSheet("color: rgb(85, 255, 255);")

    def init_quad_table(self, devs, calc_obj):
        self.ui.w_q_table.init_table(devs=devs, device_iface=DeviceUI, calc_obj=calc_obj, spin_params=[-5000, 5000, 0.1])

    #def init_s2e_table(self, sections, check_header="", phys_proc_names=[]):
    #    self.ui.w_s2e.init_table(sections=sections, check_header=check_header,
    #                             phys_proc_names=phys_proc_names)

    #def init_cavs_table(self, devs, calc_obj):
    #    self.ui.w_cav_table.init_table(devs=devs, device_iface=DeviceUI, calc_obj=calc_obj,
    #                                   spin_params_phi=[-180, 180, 0.1],
    #                                   spin_params_v=[0, 1, 0.001])


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
