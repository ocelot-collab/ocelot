# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'app/ui_forms/widgets/ebeam_table_ui.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Ebeam_Table(object):
    def setupUi(self, Ebeam_Table):
        Ebeam_Table.setObjectName("Ebeam_Table")
        Ebeam_Table.resize(700, 500)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Ebeam_Table.sizePolicy().hasHeightForWidth())
        Ebeam_Table.setSizePolicy(sizePolicy)
        self.gridLayout_2 = QtWidgets.QGridLayout(Ebeam_Table)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.table = QtWidgets.QTableWidget(Ebeam_Table)
        self.table.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.setRowCount(0)
        self.table.setColumnCount(2)
        self.table.setObjectName("table")
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.table.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        self.table.setHorizontalHeaderItem(1, item)
        self.table.verticalHeader().setVisible(False)
        self.gridLayout_2.addWidget(self.table, 0, 1, 1, 1)

        self.retranslateUi(Ebeam_Table)
        QtCore.QMetaObject.connectSlotsByName(Ebeam_Table)

    def retranslateUi(self, Ebeam_Table):
        _translate = QtCore.QCoreApplication.translate
        Ebeam_Table.setWindowTitle(_translate("Ebeam_Table", "Main parameters"))
        item = self.table.horizontalHeaderItem(0)
        item.setText(_translate("Ebeam_Table", "Parameter"))
        item = self.table.horizontalHeaderItem(1)
        item.setText(_translate("Ebeam_Table", "Value"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Ebeam_Table = QtWidgets.QWidget()
    ui = Ui_Ebeam_Table()
    ui.setupUi(Ebeam_Table)
    Ebeam_Table.show()
    sys.exit(app.exec_())

