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
        self.dockWidget = QtWidgets.QDockWidget(Widget)
        self.dockWidget.setObjectName("dockWidget")
        self.dockWidgetContents = QtWidgets.QWidget()
        self.dockWidgetContents.setObjectName("dockWidgetContents")
        self.gridLayout = QtWidgets.QGridLayout(self.dockWidgetContents)
        self.gridLayout.setContentsMargins(11, 11, 11, 11)
        self.gridLayout.setSpacing(6)
        self.gridLayout.setObjectName("gridLayout")
        self.pb_reset = QtWidgets.QPushButton(self.dockWidgetContents)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.pb_reset.setFont(font)
        self.pb_reset.setStyleSheet("color: rgb(255, 128, 0); font-size: 10pt;")
        self.pb_reset.setObjectName("pb_reset")
        self.gridLayout.addWidget(self.pb_reset, 1, 1, 1, 1)
        self.cb_realtime = QtWidgets.QCheckBox(self.dockWidgetContents)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.cb_realtime.setFont(font)
        self.cb_realtime.setObjectName("cb_realtime")
        self.gridLayout.addWidget(self.cb_realtime, 1, 0, 1, 1)
        self.pb_calculate = QtWidgets.QPushButton(self.dockWidgetContents)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.pb_calculate.setFont(font)
        self.pb_calculate.setStyleSheet("color: rgb(255, 128, 0); font-size: 10pt;")
        self.pb_calculate.setObjectName("pb_calculate")
        self.gridLayout.addWidget(self.pb_calculate, 1, 2, 1, 1)
        self.tableWidget = QtWidgets.QTableWidget(self.dockWidgetContents)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)
        self.gridLayout.addWidget(self.tableWidget, 0, 0, 1, 3)
        self.dockWidget.setWidget(self.dockWidgetContents)
        self.gridLayout_2.addWidget(self.dockWidget, 0, 0, 1, 3)

        self.retranslateUi(Widget)
        QtCore.QMetaObject.connectSlotsByName(Widget)

    def retranslateUi(self, Widget):
        _translate = QtCore.QCoreApplication.translate
        Widget.setWindowTitle(_translate("Widget", "Widget"))
        self.pb_reset.setText(_translate("Widget", "Reset"))
        self.cb_realtime.setText(_translate("Widget", "RealTime"))
        self.pb_calculate.setText(_translate("Widget", "Calculate"))

