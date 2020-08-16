# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_ocl_table.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Widget(object):
    def setupUi(self, Widget):
        Widget.setObjectName("Widget")
        Widget.resize(400, 300)
        self.gridLayout_2 = QtWidgets.QGridLayout(Widget)
        self.gridLayout_2.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_2.setSpacing(6)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.cb_realtime = QtWidgets.QCheckBox(Widget)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.cb_realtime.setFont(font)
        self.cb_realtime.setObjectName("cb_realtime")
        self.gridLayout_2.addWidget(self.cb_realtime, 1, 0, 1, 1)
        self.pb_reset = QtWidgets.QPushButton(Widget)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.pb_reset.setFont(font)
        self.pb_reset.setStyleSheet("color: rgb(255, 128, 0); font-size: 10pt;")
        self.pb_reset.setObjectName("pb_reset")
        self.gridLayout_2.addWidget(self.pb_reset, 1, 1, 1, 1)
        self.pb_calculate = QtWidgets.QPushButton(Widget)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.pb_calculate.setFont(font)
        self.pb_calculate.setStyleSheet("color: rgb(255, 128, 0); font-size: 10pt;")
        self.pb_calculate.setObjectName("pb_calculate")
        self.gridLayout_2.addWidget(self.pb_calculate, 1, 2, 1, 1)
        self.tableWidget = QtWidgets.QTableWidget(Widget)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)
        self.gridLayout_2.addWidget(self.tableWidget, 0, 0, 1, 3)

        self.retranslateUi(Widget)
        QtCore.QMetaObject.connectSlotsByName(Widget)

    def retranslateUi(self, Widget):
        _translate = QtCore.QCoreApplication.translate
        Widget.setWindowTitle(_translate("Widget", "Widget"))
        self.cb_realtime.setText(_translate("Widget", "RealTime"))
        self.pb_reset.setText(_translate("Widget", "Reset"))
        self.pb_calculate.setText(_translate("Widget", "Calculate"))

