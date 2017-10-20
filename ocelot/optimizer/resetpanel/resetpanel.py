#!/usr/local/lcls/package/python/current/bin/python
"""
PYQT interface for running OCELOT simplex optimization.
Created as a QT widget for use in other applications as well.
Tyler Cope, 2016
The class was modified and was introduced new methods.
S. Tomin, 2017
"""

import sys

from PyQt5.QtWidgets import QApplication, QFrame
from PyQt5 import QtGui, QtCore, uic
from ocelot.optimizer.mint.opt_objects import *

from ocelot.optimizer.resetpanel.UIresetpanel import Ui_Form

sys.path.append("..")


class ResetpanelWindow(QFrame):
    """
    Main GUI class for the resetpanel.
    """

    def __init__(self, parent=None):

        # initialize
        QFrame.__init__(self)
        self.ui = Ui_Form()
        self.ui.setupUi(self)

        # blank data
        self.pvs = []
        self.devices = []
        self.startValues = {}
        self.pv_objects = {}

        # button connections
        self.ui.updateReference.clicked.connect(self.updateReference)
        self.ui.resetAll.clicked.connect(self.launchPopupAll)

        # fast timer start
        self.trackTimer = QtCore.QTimer()
        self.trackTimer.timeout.connect(self.updateCurrentValues)
        self.trackTimer.start(500)  # refresh every 100 ms

        # dark theme
        self.loadStyleSheet()

    def loadStyleSheet(self):
        """ Load in the dark theme style sheet. """
        try:
            self.cssfile = "style.css"
            with open(self.cssfile, "r") as f:
                self.setStyleSheet(f.read())
        except IOError:
            print('No style sheet found!')

    # def getPvList(self, pvs_in=None):
    #    """
    #    Method to build a pv list from file.
    #
    #    Can get passed string filename to parse text into a pv list, ex.
    #            PV_1
    #            PV_2
    #            ...
    #            PV_N
    #
    #    Entrys with at '#' in the file are ignored
    #    Alternativly can be passed a list of pv strings for build the pv list
    #    Saves the PV list as a class varable
    #
    #    Args:
    #            pvs_in: Can be either List of pvs, or string filename
    #
    #    """
    #
    #    if not pvs_in:
    #        return
    #
    #    if type(pvs_in) != list:
    #        self.pvs = []
    #        for line in open(pvs_in):
    #            l = line.rstrip('\n')
    #            if l[0] == '#': #exclude commented PVs
    #                continue
    #            self.pvs.append(str(l))
    #    else:
    #        self.pvs = pvs_in
    #
    #    self.devices = self.create_devices(pvs=self.pvs)
    #
    #    self.getStartValues()
    #    self.initTable()

    def getStartValues(self):
        """ Initializes start values for the PV list. """
        for dev in self.devices:
            try:
                self.startValues[dev.eid] = dev.get_value()
            except:
                self.startValues[dev.eid] = None
                print("Get Start Value: ", dev.eid, " not working")
                # print(self.startValues[dev.eid])
                # self.pv_objects[pv].add_callback(callback=self.PvGetCallBack)

    def updateReference(self):
        """Updates reference values for all PVs on button click."""
        self.ui.updateReference.setText("Getting vals...")
        self.getStartValues()
        for row in range(len(self.pvs)):
            pv = self.pvs[row]
            self.ui.tableWidget.setItem(row, 1, QtGui.QTableWidgetItem(str(np.around(self.startValues[pv], 4))))
        self.ui.updateReference.setText("Update Reference")

    def initTable(self):
        """ Initialize the UI table object """
        headers = ["PVs", "Reference Value", "Current Value"]
        self.ui.tableWidget.setColumnCount(len(headers))
        self.ui.tableWidget.setHorizontalHeaderLabels(headers)
        self.ui.tableWidget.setEditTriggers(QtGui.QAbstractItemView.NoEditTriggers)  # No user edits on talbe
        self.ui.tableWidget.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
        for row in range(len(self.pvs)):

            self.ui.tableWidget.setRowCount(row + 1)
            pv = self.pvs[row]
            # put PV in the table
            self.ui.tableWidget.setItem(row, 0, QtGui.QTableWidgetItem(str(pv)))
            #self.ui.tableWidget.item(row, 0).setTextColor(QtGui.QColor(0, 255, 255))
            self.ui.tableWidget.item(row, 0).setForeground(QtGui.QColor(0, 255, 255))
            tip = "/".join(str(pv).split("/")[-2:])
            self.ui.tableWidget.item(row, 0).setToolTip(tip)
            # self.ui.tableWidget.item(row, 0).setFont(font)
            # put start val in
            s_val = self.startValues[pv]
            if s_val != None:
                s_val = np.around(s_val, 4)
            self.ui.tableWidget.setItem(row, 1, QtGui.QTableWidgetItem(str(s_val)))

            # change font size
            # font = QtGui.QFont()
            # font.setPointSize(12)
            # self.ui.tableWidget.item(row, 1).setFont(font)

            # self.pv_objects[pv].run_callbacks()#initialize in the pvs current value
            # print("init", self.ui.tableWidget.item(row, 0))

    # update the table on PV change callback
    #
    # REMOVED BECAUSE CALLBACK CAUSED SEG FAUTLS.
    #
    # def PvGetCallBack(self,**kw):

    #        #get pv info from the callback kw arg
    #        val=kw['value']
    #        pv=kw['pvname']
    #
    #        #set current string to red at 0.5 percent difference from initial value
    #        percent = 0.005

    #
    #        #find the difference from PV start to current, decide to change table color
    #        try:
    #                #get row for callback PV in the gui table, set string to pv value
    #                row = self.pvs.index(str(pv))
    #                self.ui.tableWidget.setItem(row,2,QtGui.QTableWidgetItem(str(val)))

    #                tol  = abs(self.startValues[pv]*percent)
    #                diff = abs(abs(self.startValues[pv]) - abs(val))
    #                if diff > tol:
    #                        self.ui.tableWidget.item(row,2).setForeground(QtGui.QColor(255,0,0))
    #                else:
    #                        self.ui.tableWidget.item(row,2).setForeground(QtGui.QColor(255,255,255))
    #                QApplication.processEvents()
    #        except:
    #                pass

    def updateCurrentValues(self):

        """
        Method to update the table on every clock cycle.
        Loops through the pv list and gets new data, then updates the Current Value column.
        Hard coded to turn Current Value column red at 0.1% differenct from Ref Value.
        It would be better to update the table on a callback, but PyEpics crashes with cb funcitons.
        """
        percent = 0.001
        self.currentValues = {}
        for row, dev in enumerate(self.devices):
            try:
                value = dev.get_value()
            except:
                # print("ERROR getting value. Device:", dev.eid)
                value = None

            if self.startValues[dev.eid] == None and value != None:
                self.startValues[dev.eid] = value

            if self.startValues[dev.eid] == None or value == None:
                self.ui.tableWidget.item(row, 5).setFlags(QtCore.Qt.NoItemFlags)
                for col in [0, 5]:
                    self.ui.tableWidget.item(row, col).setBackground(QtGui.QColor(255, 0, 0))  # red

                if self.startValues[dev.eid] == None:
                    self.ui.tableWidget.setItem(row, 1, QtGui.QTableWidgetItem(str("None")))
                    self.ui.tableWidget.item(row, 1).setBackground(QtGui.QColor(255, 0, 0))  # red
                else:
                    self.ui.tableWidget.setItem(row, 1, QtGui.QTableWidgetItem(str(np.around(value, 4))))
                    self.ui.tableWidget.item(row, 1).setBackground(QtGui.QColor(89, 89, 89))  # grey

                if value == None:
                    self.ui.tableWidget.setItem(row, 2, QtGui.QTableWidgetItem(str("None")))
                    self.ui.tableWidget.item(row, 2).setBackground(QtGui.QColor(255, 0, 0))  # red
                else:
                    self.ui.tableWidget.setItem(row, 2,
                                                QtGui.QTableWidgetItem(str(np.around(self.currentValues[dev.eid], 4))))
                    self.ui.tableWidget.item(row, 2).setBackground(QtGui.QColor(89, 89, 89))  # grey

                continue

            # if value out of the limits
            if dev.check_limits(value):
                for col in [3, 4]:
                    spin_box = self.ui.tableWidget.cellWidget(row, col)
                    spin_box.setStyleSheet("color: yellow; font-size: 16px; background-color:red;")

            else:
                for col in [3, 4]:
                    spin_box = self.ui.tableWidget.cellWidget(row, col)
                    if col == 3:
                        spin_box.setStyleSheet("color: rgb(153,204,255); font-size: 16px; background-color:#595959;")
                    if col == 4:
                        spin_box.setStyleSheet("color: rgb(255,0,255); font-size: 16px; background-color:#595959;")

            pv = dev.eid

            self.currentValues[pv] = value  # dev.get_value()
            self.ui.tableWidget.setItem(row, 2, QtGui.QTableWidgetItem(str(np.around(self.currentValues[pv], 4))))
            # print(self.currentValues[pv])
            tol = abs(self.startValues[pv] * percent)
            diff = abs(abs(self.startValues[pv]) - abs(self.currentValues[pv]))
            if diff > tol:
                self.ui.tableWidget.item(row, 2).setForeground(QtGui.QColor(255, 101, 101))  # red
            else:
                self.ui.tableWidget.item(row, 2).setForeground(QtGui.QColor(255, 255, 255))  # white

            self.ui.tableWidget.item(row, 5).setFlags(
                QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)

            for col in [0, 1, 2, 5]:
                # self.ui.tableWidget.item(row, col).setFlags(
                #    QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsUserCheckable | QtCore.Qt.ItemIsEnabled)
                self.ui.tableWidget.item(row, col).setBackground(QtGui.QColor(89, 89, 89))

        QApplication.processEvents()

    def resetAll(self):
        """Set all PVs back to their reference values."""
        for dev in self.devices:
            val = self.startValues[dev.eid]
            dev.set_value(val)  # epics.caput(pv,val)

    def launchPopupAll(self):
        """Launches the ARE YOU SURE popup window for pv reset."""
        self.ui_check = uic.loadUi("UIareyousure.ui")
        self.ui_check.exit.clicked.connect(self.ui_check.close)
        self.ui_check.reset.clicked.connect(self.resetAll)
        self.ui_check.reset.clicked.connect(self.ui_check.close)
        self.ui_check.show()


def main():
    """
    Main functino to open a resetpanel GUI.
    If passed a file name, will try and load PV list from that file.
    Otherwise defaults to a file in the base directory with pre-loaded common tuned PVs.
    """
    try:  # try to get a pv list file name from commandline arg
        pvs = sys.argv[1]
    except:
        pvs = "./lclsparams"

    app = QApplication(sys.argv)
    window = ResetpanelWindow()
    window.setWindowIcon(QtGui.QIcon('/usr/local/lcls/tools/python/toolbox/py_logo.png'))
    window.getPvList(pvs)
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

