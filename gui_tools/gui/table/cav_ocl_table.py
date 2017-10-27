import sys
from PyQt5.QtWidgets import QCheckBox, QHBoxLayout, QHeaderView, QApplication,QMenu, QWidget, QAction, QTableWidget, QTableWidgetItem, QDoubleSpinBox

from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from gui.table.ui_ocl_table import *
from gui.table.v_ocl_table import *
import numpy as np


class CavOclTable(VOclTable):
    def __init__(self, parent=None):
        super().__init__()
        self.ui.cb_realtime.stateChanged.connect(self.realtime)
        self.ui.cb_realtime.setChecked(False)
        self.ui.pb_reset.clicked.connect(self.reset)


    def init_table(self, devs, device_iface=None, calc_obj=None, spin_params_phi=[-180, 180, 0.1],
                   spin_params_v=[0, 1, 0.001], check_box=False):
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

        h_headers = ["ID", "Phi.ref.", "Phi.cur.", 'V.ref.', "V.cur."]
        if check_box:
            h_headers += "status"
        self.ui.tableWidget.setColumnCount(len(h_headers))
        self.ui.tableWidget.setHorizontalHeaderLabels(h_headers)

        self.ui.tableWidget.setRowCount(0)
        for row, dev in enumerate(devs):
            #print(dev)
            self.ui.tableWidget.setRowCount(row + 1)
            pv = dev.id
            # put ID in the table
            self.ui.tableWidget.setItem(row, 0, QTableWidgetItem(str(pv)))

            # put init val in
            val_phi = np.round(dev.phi, 4)
            self.ui.tableWidget.setItem(row, 1, QTableWidgetItem(str(val_phi)))

            spin_box_phi = self.create_spin_box(spin_params_phi)
            spin_box_phi.setValue(dev.phi)

            val_v = np.round(dev.v, 4)
            self.ui.tableWidget.setItem(row, 3, QTableWidgetItem(str(val_v)))
            spin_box_v = self.create_spin_box(spin_params_v)
            spin_box_v.setValue(dev.v)

            if calc_obj != None:
                spin_box_phi.valueChanged.connect(calc_obj)
                spin_box_v.valueChanged.connect(calc_obj)
            self.ui.tableWidget.setCellWidget(row, 2, spin_box_phi)
            self.ui.tableWidget.setCellWidget(row, 4, spin_box_v)

            if check_box:
                checkBoxItem = QTableWidgetItem()
                # checkBoxItem.setBackgroundColor(QtGui.QColor(100,100,150))
                checkBoxItem.setCheckState(QtCore.Qt.Checked)
                flags = checkBoxItem.flags()
                # print("FLAG", flags)
                # flags != flags
                checkBoxItem.setFlags(flags)
                self.ui.tableWidget.setItem(row, 5, checkBoxItem)
                dev.row = row

            if device_iface != None:
                ui = device_iface()
                ui.tableWidget = self.ui.tableWidget
                ui.row = row
                ui.col = 2
                dev.ui = ui

        self.ui.tableWidget.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.realtime()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = VOclTable()
    window.show()
    sys.exit(app.exec_())
