#!/usr/local/lcls/package/python/current/bin/python
"""
Subclass of ResetpanelWindow, used to make a version with check box selction

Tyler Cope, 2016

The class was modified and was introduced new methods.

S. Tomin, 2017
"""

from ocelot.optimizer.resetpanel.resetpanel import ResetpanelWindow
from PyQt5.QtWidgets import QApplication, QFrame, QPushButton, QTableWidget
from PyQt5 import QtGui, QtCore, Qt, uic
from PyQt5.QtGui import QClipboard
import sys
import time
#from ocelot.optimizer.mint.opt_objects import *
from ocelot.optimizer.mint import opt_objects as obj
from ocelot.optimizer.mint.lcls_interface import *


class customTW(QTableWidget):
    """
    Subclass the tablewidget to add a custom PV on middle click.
    """
    def __init__(self, parent):
        """Init for the table widget, nothing special here"""
        QTableWidget.__init__(self, parent)
        self.parent = parent

    def mouseReleaseEvent(self, evt):
        """
        Grabs string from the clipboard if middle button click.

        Tries to then add the PV to the parent GUI pv lists.
        Calls parent.addPv() method.

        Args:
                evt (QEvent): Event object passed in from QT
        """
        button = evt.button()
        #print("TW Release event: ", button, self.parent.enableMiddleClick)
        if (button == 4) and (self.parent.enableMiddleClick):

            pv = QtGui.QApplication.clipboard().text()
            self.parent.addPv(pv)
            #print(QtGui.QApplication.clipboard().text())
        else:
            QTableWidget.mouseReleaseEvent(self, evt)

    def mouseReleaseEvent_mclick(self, evt):
        """
        Grabs string from the clipboard if middle button click.

        Tries to then add the PV to the parent GUI pv lists.
        Calls parent.addPv() method.

        Args:
                evt (QEvent): Event object passed in from QT
        """
        button = evt.button()
        if (button == 4) and (self.parent.enableMiddleClick):
            pv = QtGui.QApplication.clipboard().text(mode=QClipboard.Selection)
            self.parent.addPv(pv)
        else:
            QTableWidget.mouseReleaseEvent(self, evt)

    def contextMenuEvent(self, event):
        if self.selectionModel().selection().indexes():
            rows = []
            for i in self.selectionModel().selection().indexes():
                row, column = i.row(), i.column()
                rows.append(row)
            self.menu = QtGui.QMenu(self)
            deleteAction = QtGui.QAction('Delete', self)
            deleteAction.triggered.connect(lambda: self.deleteSlot(rows))

            self.menu.addAction(deleteAction)

            #editAction = QtGui.QAction('Edit', self)
            #self.menu.addAction(editAction)
            # add other required actions
            self.menu.popup(QtGui.QCursor.pos())


    def deleteSlot(self, rows):
        print ("delete rows called", rows)
        for row in rows[::-1]:
            self.parent.ui.tableWidget.removeRow(row)

        table = self.parent.get_state()
        pvs = table["id"]
        self.parent.getPvList(pvs)
        self.parent.uncheckBoxes()
        self.parent.set_state(table)


class ResetpanelBoxWindow(ResetpanelWindow):
    """
    The main GUI class to add in checkboxes, subclassed from the ResetPanelWindow class.
    """
    def __init__(self,  parent=None):
        """
        Init method, adds in some new UI changes.

        Adds two buttons, to check and uncheck selected rows.
        Add in subclassed table to enable middle click PV add.
        """

        #initialize

        ResetpanelWindow.__init__(self)

        self.check = QPushButton(self)
        self.check.setText('Check')
        self.uncheck = QPushButton(self)
        self.uncheck.setText('Uncheck')
        self.ui.horizontalLayout.addWidget(self.check)
        self.ui.horizontalLayout.addWidget(self.uncheck)
        self.check.clicked.connect(lambda: self.getRows(2))
        self.uncheck.clicked.connect(lambda: self.getRows(0))
        self.mi = None
        self.dp = None

        #make button text bigger
        #self.check.setStyleSheet('font-size: 18pt; font-family: Courier;')

        #enable middle click method
        self.enableMiddleClick = True

        #make the custom table for middle click
        self.ui.tableWidget.setParent(None) #remove old table
        self.ui.tableWidget = customTW(self) # make new widget
        self.ui.gridLayout.addWidget(self.ui.tableWidget,0,0)
        #self.ui.tableWidget.itemClicked.connect(self.con)


    def mouseReleaseEvent(self, evt):
        """
        Get PV coppied from the system clipboard on button release

        Button 4 is the middle click button.
        """
        button = evt.button()
        #print("Release event: ", button, self.parent.enableMiddleClick)
        if (button == 4) and (self.enableMiddleClick):
            pv = QtGui.QApplication.clipboard().text(mode=QClipboard.Selection)
            self.addPv(pv)

    def set_parent(self, parent):
        self.parent = parent

    def addPv(self, pv):
        """
        Add another PV to the GUI on middle click.

        Args:
                pv (str): String name of the PV to add
        """
        pv = str(pv)
        if pv in self.pvs:
            print ("PV already in list")
            return
        try:
            #print("try to create dev")
            dev = self.parent.create_devices(pvs=[pv])[0]#obj.TestDevice(eid=pv)

        except:
            print ("bad string")
            return

        state = dev.state()
        print("state=", state)
        if state:
            self.pvs.append(pv)
            self.devices.append(dev)
            #self.pv_objects[pv].add_callback(callback=self.PvGetCallBack)
            self.getStartValues()
        table = self.get_state()
        self.initTable()
        self.addCheckBoxes()
        self.uncheckBoxes()
        self.set_state(table)



    def get_devices(self, pvs):
        d_pvs = [dev.eid for dev in self.devices]
        inxs = [d_pvs.index(pv) for pv in pvs]
        return [self.devices[inx] for inx in inxs]

    def set_machine_interface(self, mi, dp):
        self.mi = mi
        self.dp = dp
        self.getPvList(self.pvs)

    def getPvList(self, pvs_in=None):
        """
        Redefine method to add in checkboxes when getting the PV list.

        Copied from resetpanel.py
        """
        print ("LOADING:")
        if pvs_in == None:
            #print ('Exiting', pvs_in)
            return

        if type(pvs_in) != list:
            self.pvs = []
            #print(pvs_in)
            for line in open(pvs_in):
                l = line.rstrip('\n')
                if l[0] == '#':
                    continue
                self.pvs.append(str(l))
        else:
            self.pvs = pvs_in

        #print ("PVS LOADED", self.pvs)
        self.devices = self.parent.create_devices(self.pvs)
        self.getStartValues()
        self.initTable()
        self.addCheckBoxes()

    def get_state(self):
        devs = {"id":[], "lims": []}
        for row in range(self.ui.tableWidget.rowCount()):
            name = str(self.ui.tableWidget.item(row, 0).text())
            devs["id"].append(name)
            devs["lims"].append([self.ui.tableWidget.cellWidget(row, 3).value(), self.ui.tableWidget.cellWidget(row, 4).value()])
        return devs

    def get_limits(self, pv):
        for row in range(self.ui.tableWidget.rowCount()):
            if pv == str(self.ui.tableWidget.item(row, 0).text()):
                lims = [self.ui.tableWidget.cellWidget(row, 3).value(), self.ui.tableWidget.cellWidget(row, 4).value()]
                return lims
        return None

    def set_state(self, table):
        for row in range(self.ui.tableWidget.rowCount()):
            pv = str(self.ui.tableWidget.item(row, 0).text())
            if pv in table["id"]:
                indx = table["id"].index(pv)
                self.ui.tableWidget.cellWidget(row, 3).setValue(table["lims"][indx][0])
                self.ui.tableWidget.cellWidget(row, 4).setValue(table["lims"][indx][1])


        #for row, dev in enumerate(table["id"]):
        #    #put PV in the table
        #    self.ui.tableWidget.item(row, 0).setText(dev)
        #    self.ui.tableWidget.cellWidget(row, 3).setValue(table["lims"][row][0])
        #    self.ui.tableWidget.cellWidget(row, 4).setValue(table["lims"][row][1])

    def getRows(self, state):
        """
        Method to set the UI checkbox state from slected rows.

        Loops though the rows and gets the selected state from the 'Active" column.
        If highlighted, check box is set the the 'state' input arg.

        Args:
                state (bool): Bool of whether the boxes should be checked or unchecked.
        """
        rows=[]
        for idx in self.ui.tableWidget.selectedIndexes():
            rows.append(idx.row())
            #item = self.ui.tableWidget.cellWidget(idx.row(), 5)
            item = self.ui.tableWidget.item(idx.row(), 5)
            if item.flags() == QtCore.Qt.NoItemFlags:
                print("item disabled")
                continue
            item.setCheckState(state)
            #print item.text()
        #print rows

    def addCheckBoxes(self):
        """
        Creats additional column in UI table for check box.

        Must be called again if the user adds another PV with middle click function.
        """
        headers = ["PVs", "Saved Val.", "Current Val.", "Min", "Max", "Active"]
        self.ui.tableWidget.setColumnCount(len(headers))
        self.ui.tableWidget.setHorizontalHeaderLabels(headers)
        header = self.ui.tableWidget.horizontalHeader()
        header.setResizeMode(0, QtGui.QHeaderView.Stretch)
        header.setResizeMode(1, QtGui.QHeaderView.ResizeToContents)
        header.setResizeMode(2, QtGui.QHeaderView.ResizeToContents)
        header.setResizeMode(3, QtGui.QHeaderView.ResizeToContents)
        header.setResizeMode(4, QtGui.QHeaderView.ResizeToContents)
        #header.setResizeMode(5, QtGui.QHeaderView.ResizeToContents)
        header.setResizeMode(5, QtGui.QHeaderView.Fixed)
        #self.ui.tableWidget.horizontalHeader().resizeSection(5, 80)
        #self.ui.tableWidget.horizontalHeader().resizeSection(3, 80)
        for row in range(len(self.pvs)+1):
            eng = QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates)
            #spin_box = QtGui.QSpinBox()
            for i in range(2):
                spin_box = QtGui.QDoubleSpinBox()
                if i == 0:
                    spin_box.setStyleSheet("color: rgb(153,204,255); font-size: 16px; background-color:#595959;")
                else:
                    spin_box.setStyleSheet("color: rgb(255,0,255); font-size: 16px; background-color:#595959;")
                spin_box.setLocale(eng)
                spin_box.setDecimals(3)
                spin_box.setMaximum(999)
                spin_box.setMinimum(-999)
                spin_box.setSingleStep(0.1)
                spin_box.setAccelerated(True)
                #spin_box.setFixedWidth(50)
                self.ui.tableWidget.setCellWidget(row, 3+i, spin_box)
                self.ui.tableWidget.resizeColumnsToContents()
                #self.ui.tableWidget.setItem(row, 3 + i, spin_box)

            #checkBoxItem = QtGui.QCheckBox()
            #checkBoxItem.setStyleSheet("background-color:#595959;")
            #self.ui.tableWidget.setCellWidget(row, 5, checkBoxItem)

            #spin_box
            checkBoxItem = QtGui.QTableWidgetItem()
            #checkBoxItem.setBackgroundColor(QtGui.QColor(100,100,150))
            checkBoxItem.setCheckState(QtCore.Qt.Checked)
            flags = checkBoxItem.flags()
            #print("FLAG", flags)
            #flags != flags
            checkBoxItem.setFlags(flags)
            self.ui.tableWidget.setItem(row, 5, checkBoxItem)
            #self.ui.tableWidget.setItem(row, 4, spin_box)


    def uncheckBoxes(self):
        """ Method to unchecked all active boxes """
        for row in range(len(self.pvs)):
            item=self.ui.tableWidget.item(row, 5)
            #item = self.ui.tableWidget.cellWidget(row, 5)
            item.setCheckState(False)

    def resetAll(self):
        """
        Resets all PVs with a box checked.

        Rewrote this function to only change selected rows, not all rows.
        """
        for row, dev in enumerate(self.devices):
            val = self.startValues[dev.eid]
            state = self.ui.tableWidget.item(row, 5).checkState()
            if state == 2:
                dev.set_value(val)
                #epics.caput(pv, val)

    def updateReference(self):
        """
        Update reference values for selected rows.

        Rewrote this function to only update slected rows, not all.
        """
        self.ui.updateReference.setText("Getting vals...")
        for row, dev in enumerate(self.devices):
            pv = self.pvs[row]
            state = self.ui.tableWidget.item(row, 5).checkState()
            print("STATE")
            if state == 2:
                self.startValues[pv] = dev.get_value()
                self.ui.tableWidget.setItem(row, 1, QtGui.QTableWidgetItem(str(np.around(self.startValues[pv], 4))))
                #self.pv_objects[pv] = epics.PV(pv)
                #self.pv_objects[pv].add_callback(callback=self.PvGetCallBack)
        self.ui.updateReference.setText("Update Reference")


    def getPvsFromCbState(self):

        """
        Gets list of all pvs that have checked boxes.

        Returns:
                List of PV strings
        """
        pvs = []
        for row in range(len(self.pvs)):
            state = self.ui.tableWidget.item(row, 5).checkState()
            if state == 2:
                pvs.append(str(self.ui.tableWidget.item(row, 0).text()))
        return pvs


    #switch string from to SLECTED
    def launchPopupAll(self):
        """
        Launches the ARE YOU SURE popup for pv reset value function.

        Rewrote to change the warning string to say "checkbox selected" instead of "all" to avoid confusion with number of devices being reset.
        """
        #self.ui_check = uic.loadUi("/home/physics/tcope/python/tools/resetpanel/UIareyousure.ui")
        self.ui_check = uic.loadUi("../optimizer/resetpanel/UIareyousure.ui")
        self.ui_check.exit.clicked.connect(self.ui_check.close)
        self.ui_check.reset.clicked.connect(self.resetAll)
        self.ui_check.reset.clicked.connect(self.ui_check.close)
        self.ui_check.label.setText("Are you sure you want to implement \nchanges to checkbox selected PVs?")
        self.ui_check.show()

def main():
    """
    Start up the main program if launch from comamnd line.
    """
    try:
        pvs = sys.argv[1]
    except:
        pvs = "./lclsparams"

    app = QApplication(sys.argv)
    window = ResetpanelBoxWindow()
    window.setWindowIcon(QtGui.QIcon('/usr/local/lcls/tools/python/toolbox/py_logo.png'))
    window.getPvList(pvs)
    window.uncheckBoxes()
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
