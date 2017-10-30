"""
S.Tomin. Interface between device obgect and GUI
"""
from PyQt5 import  QtCore, QtGui
import numpy as np

class SectionUI:
    def __init__(self, ui=None):
        self.tableWidget = None
        self.row = 0
        self.col = 0
        self.alarm = False
        self.v_headers = []
        # self.tableWidget.item(self.row, 1).setBackground(QtGui.QColor(89, 89, 89))  # grey
        # self.tableWidget.item(self.row, 1).setBackground(QtGui.QColor(255, 0, 0))  # red
        # self.tableWidget.item(self.row, 1).setText(str(x))

    def get_status(self, proc=""):
        for row, vheader in enumerate(self.v_headers):
            if proc == vheader:
                return self.state_row(row)
        print("SectionUI: get_status-> No proc")
        return False

    def uncheck_proc(self, proc=""):
        for row, vheader in enumerate(self.v_headers):
            if proc == vheader:
                self.uncheck(row)
                return True
        print("SectionUI: uncheck_proc-> No proc")
        return False

    def check_proc(self, proc=""):
        for row, vheader in enumerate(self.v_headers):
            if proc == vheader:
                self.check(row)
                return True
        print("SectionUI: uncheck_proc-> No proc")
        return False

    def disable_proc(self, proc=""):
        for row, vheader in enumerate(self.v_headers):
            if proc == vheader:
                self.disable(row)
                return True
        print("SectionUI: disable_proc-> No proc")
        return False

    def enable_proc(self, proc=""):
        for row, vheader in enumerate(self.v_headers):
            if proc == vheader:
                self.enable(row)
                return True
        print("SectionUI: enable_proc-> No proc")
        return False

    def uncheck(self, row):
        item = self.tableWidget.item(row, self.col)
        item.setCheckState(False)

    def check(self, row):
        item = self.tableWidget.item(row, self.col)
        item.setCheckState(QtCore.Qt.Checked)

    def enable(self, row):
        checkBoxItem = self.tableWidget.item(row, self.col)
        checkBoxItem.setFlags(QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsUserCheckable)

    def disable(self, row):
        checkBoxItem = self.tableWidget.item(row, self.col)
        checkBoxItem.setCheckState(False)
        checkBoxItem.setFlags(QtCore.Qt.ItemIsSelectable)

    def calculated(self, flag):
        checkBoxItem = self.tableWidget.item(0, self.col)
        if flag:
            checkBoxItem.setBackground(QtGui.QColor(0, 255, 0))
        else:
            checkBoxItem.setBackground(QtGui.QColor(100, 100, 150))

    #def get_status(self, row):
    #    item = self.tableWidget.item(row, self.col)
    #    if item.checkState() == QtCore.Qt.Checked:
    #        return True
    #    else:
    #        return False

    def state_row(self, row):
        item = self.tableWidget.item(row, self.col)
        state = item.checkState()
        return state

    def set_hide_row(self, hide):
        #if hide:
        #    self.uncheck()
        #else:
        #    self.check()
        self.tableWidget.setRowHidden(self.row, hide)


class DeviceUI:
    def __init__(self, ui=None):
        self.tableWidget = None
        self.row = 0
        self.col = 0
        self.alarm = False

    def get_value(self):
        return self.tableWidget.cellWidget(self.row, self.col).value()

    def set_value(self, val):
        self.tableWidget.cellWidget(self.row, self.col).setValue(val)

    def set_init_value(self, val):
        val = np.round(val, 4) # "{:1.4e}".format(val)
        self.tableWidget.item(self.row, 1).setText(str(val))

    def get_init_value(self):
        return float(self.tableWidget.item(self.row, 1).text())

    def uncheck(self):
        item = self.tableWidget.item(self.row, 3)
        item.setCheckState(False)

    def check(self):
        item = self.tableWidget.item(self.row, 3)
        item.setCheckState(QtCore.Qt.Checked)

    def state(self):
        item = self.tableWidget.item(self.row, 3)
        state = item.checkState()
        return state

    def set_alarm(self, flag):
        if flag:
            self.tableWidget.item(self.row, 0).setBackground(QtGui.QColor(255, 0, 0))  # red
        else:
            self.tableWidget.item(self.row, 0).setBackground(QtGui.QColor(89, 89, 89))  # grey

    def check_values(self, val, lims, warn=False):
        if warn:
            self.tableWidget.item(self.row, 0).setBackground(QtGui.QColor(255, 255, 0))  # yellow
        else:
            #print("grey")
            self.tableWidget.item(self.row, 0).setBackground(QtGui.QColor(89, 89, 89))  # grey
        self.alarm = False
        if not(lims[0]<= val <= lims[1]):
            self.tableWidget.item(self.row, 0).setBackground(QtGui.QColor(255, 0, 0))  # red
            self.alarm = True

    def value_was_changed(self, flag):
        if flag:
            # self.r_items[elem.ui.row].setBrush(pg.mkBrush("r"))
            self.tableWidget.item(self.row, 1).setForeground(QtGui.QColor(255, 101, 101))  # red
        else:
            self.tableWidget.item(self.row, 1).setForeground(QtGui.QColor(255, 255, 255))  # white

    def set_hide(self, hide):
        #if hide and uncheck:
        #    self.uncheck()
        #else:
        #    self.check()
        self.tableWidget.setRowHidden(self.row, hide)
        
class CavityUI:
    def __init__(self, ui=None):
        self.tableWidget = None
        self.row = 0
        self.col = 0
        self.alarm = False

    def get_volt(self):
        return self.tableWidget.cellWidget(self.row, 4).value()

    def get_phi(self):
        return self.tableWidget.cellWidget(self.row, 2).value()

    def set_volt(self, val):
        self.tableWidget.cellWidget(self.row, 4).setValue(val)

    def set_phi(self, val):
        self.tableWidget.cellWidget(self.row, 2).setValue(val)

    def set_init_volt(self, val):
        val = np.round(val, 4) # "{:1.4e}".format(val)
        self.tableWidget.item(self.row, 3).setText(str(val))
        

    def set_init_phi(self, val):
        val = np.round(val, 4) # "{:1.4e}".format(val)
        self.tableWidget.item(self.row, 1).setText(str(val))

    def get_init_phi(self):
        return float(self.tableWidget.item(self.row, 1).text())

    def get_init_volt(self):
        return float(self.tableWidget.item(self.row, 3).text())

    def uncheck(self):
        item = self.tableWidget.item(self.row, 3)
        item.setCheckState(False)

    def check(self):
        item = self.tableWidget.item(self.row, 3)
        item.setCheckState(QtCore.Qt.Checked)

    def state(self):
        item = self.tableWidget.item(self.row, 3)
        state = item.checkState()
        return state

    def set_alarm(self, flag):
        if flag:
            self.tableWidget.item(self.row, 0).setBackground(QtGui.QColor(255, 0, 0))  # red
        else:
            self.tableWidget.item(self.row, 0).setBackground(QtGui.QColor(89, 89, 89))  # grey

    def check_values(self, val, lims, warn=False):
        if warn:
            self.tableWidget.item(self.row, 0).setBackground(QtGui.QColor(255, 255, 0))  # yellow
        else:
            #print("grey")
            self.tableWidget.item(self.row, 0).setBackground(QtGui.QColor(89, 89, 89))  # grey
        self.alarm = False
        if not(lims[0]<= val <= lims[1]):
            self.tableWidget.item(self.row, 0).setBackground(QtGui.QColor(255, 0, 0))  # red
            self.alarm = True

    def phi_was_changed(self, flag):
        if flag:
            # self.r_items[elem.ui.row].setBrush(pg.mkBrush("r"))
            self.tableWidget.item(self.row, 1).setForeground(QtGui.QColor(255, 101, 101))  # red
        else:
            self.tableWidget.item(self.row, 1).setForeground(QtGui.QColor(255, 255, 255))  # white

    def volt_was_changed(self, flag):
        if flag:
            # self.r_items[elem.ui.row].setBrush(pg.mkBrush("r"))
            self.tableWidget.item(self.row, 3).setForeground(QtGui.QColor(255, 101, 101))  # red
        else:
            self.tableWidget.item(self.row, 3).setForeground(QtGui.QColor(255, 255, 255))  # white

    def set_hide(self, hide):
        #if hide and uncheck:
        #    self.uncheck()
        #else:
        #    self.check()
        self.tableWidget.setRowHidden(self.row, hide)