import sys
from PyQt5.QtWidgets import QCheckBox, QHBoxLayout, QHeaderView, QApplication,QMenu, QWidget, QAction, QTableWidget, QTableWidgetItem, QDoubleSpinBox

from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from gui.table.ui_ocl_table import *
from gui.table.ocl_table import *
import numpy as np


class VOclTable(OclTable):
    def __init__(self, parent=None):
        super().__init__()
        self.ui.cb_realtime.stateChanged.connect(self.realtime)
        self.ui.cb_realtime.setChecked(True)
        self.ui.pb_reset.clicked.connect(self.reset)

    def block_spin_box(self, flag=False):
        cols = self.ui.tableWidget.columnCount()
        for row in range(self.ui.tableWidget.rowCount()):
            self.ui.tableWidget.cellWidget(row, 2).blockSignals(flag)
            if cols > 4 and self.ui.tableWidget.cellWidget(row, 4).__class__ == QDoubleSpinBox:
                self.ui.tableWidget.cellWidget(row, 4).blockSignals(flag)

    def realtime(self):
        if self.ui.cb_realtime.checkState():
            self.block_spin_box(flag=False)
            #self.ui.tableWidget.blockSignals(True)
        else:
            self.block_spin_box(flag=True)

    def reset(self):
        self.block_spin_box(flag=True)
        for row in range(self.ui.tableWidget.rowCount()):
            val = float(self.ui.tableWidget.item(row, 1).text())
            self.ui.tableWidget.cellWidget(row, 2).setValue(val)
        self.realtime()
        self.ui.pb_calculate.click()



    #def value_was_changed(self, flag):
    #    if flag:
    #        # self.r_items[elem.ui.row].setBrush(pg.mkBrush("r"))
    #        self.ui.tableWidget.item(self.row, 1).setForeground(QtGui.QColor(255, 101, 101))  # red
    #    else:
    #        self.ui.tableWidget.item(self.row, 1).setForeground(QtGui.QColor(255, 255, 255))  # white

    def create_spin_box(self, spin_params=[-5000, 5000, 5]):
        eng = QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates)

        spin_box = QDoubleSpinBox()
        spin_box.setStyleSheet("color: #b1b1b1; font-size: 16px; background-color:#595959; border: 2px solid #b1b1b1")
        spin_box.setLocale(eng)
        spin_box.setDecimals(4)
        spin_box.setMaximum(spin_params[1])
        spin_box.setMinimum(spin_params[0])
        spin_box.setSingleStep(spin_params[2])
        spin_box.setAccelerated(True)
        return spin_box

    def init_table(self, devs, dev_param="k1", device_iface=None, calc_obj=None, spin_params=[-5000, 5000, 5], check_box=False):
        """
        Initialize the UI table object

        :param devs:
        :param device_iface:
        :param calc_obj:
        :param spin_params:
        :param check_box:
        :return:
        """
        self.calc_obj = calc_obj

        h_headers = ["ID", "ref.val.", "cur.val."]
        if check_box:
            h_headers += "status"
        self.ui.tableWidget.setColumnCount(len(h_headers))
        self.ui.tableWidget.setHorizontalHeaderLabels(h_headers)

        self.spin_boxes = []
        self.ui.tableWidget.setRowCount(0)
        for row, dev in enumerate(devs):
            #print(dev)
            eng = QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates)
            self.ui.tableWidget.setRowCount(row + 1)
            pv = dev.id
            # put PV in the table
            self.ui.tableWidget.setItem(row, 0, QTableWidgetItem(str(pv)))
            # put start val in
            if dev_param not in dev.__dict__.keys():
                print("Parameter ", dev_param, " is not in the device!")
                continue
            val = np.round(dev.__dict__[dev_param], 4)
            self.ui.tableWidget.setItem(row, 1, QTableWidgetItem(str(val)))

            spin_box = self.create_spin_box(spin_params=spin_params)
            spin_box.setValue(dev.__dict__[dev_param])

            if calc_obj != None:
                spin_box.valueChanged.connect(calc_obj)
            self.ui.tableWidget.setCellWidget(row, 2, spin_box)
            self.spin_boxes.append(spin_box)

            if check_box:
                checkBoxItem = QTableWidgetItem()
                # checkBoxItem.setBackgroundColor(QtGui.QColor(100,100,150))
                checkBoxItem.setCheckState(QtCore.Qt.Checked)
                flags = checkBoxItem.flags()
                # print("FLAG", flags)
                # flags != flags
                checkBoxItem.setFlags(flags)
                self.ui.tableWidget.setItem(row, 3, checkBoxItem)

                dev.row = row
            if device_iface != None:
                ui = device_iface()
                ui.tableWidget = self.ui.tableWidget
                ui.row = row
                ui.col = 2
                dev.ui = ui

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = VOclTable()
    window.show()
    sys.exit(app.exec_())
