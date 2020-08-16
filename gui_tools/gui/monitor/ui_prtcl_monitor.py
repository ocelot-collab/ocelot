# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ui_prtcl_monitor.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Widget(object):
    def setupUi(self, Widget):
        Widget.setObjectName("Widget")
        Widget.resize(442, 372)
        self.gridLayout_2 = QtWidgets.QGridLayout(Widget)
        self.gridLayout_2.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_2.setSpacing(6)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.groupBox_2 = QtWidgets.QGroupBox(Widget)
        self.groupBox_2.setObjectName("groupBox_2")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.groupBox_2)
        self.gridLayout_3.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_3.setSpacing(6)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.cb_slice_params = QtWidgets.QComboBox(self.groupBox_2)
        self.cb_slice_params.setObjectName("cb_slice_params")
        self.gridLayout_3.addWidget(self.cb_slice_params, 0, 0, 1, 1)
        self.gridLayout_2.addWidget(self.groupBox_2, 1, 0, 1, 1)
        self.w_monitor = QtWidgets.QWidget(Widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.w_monitor.sizePolicy().hasHeightForWidth())
        self.w_monitor.setSizePolicy(sizePolicy)
        self.w_monitor.setObjectName("w_monitor")
        self.gridLayout_2.addWidget(self.w_monitor, 0, 0, 1, 4)
        self.groupBox = QtWidgets.QGroupBox(Widget)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout.setContentsMargins(11, 11, 11, 11)
        self.gridLayout.setSpacing(6)
        self.gridLayout.setObjectName("gridLayout")
        self.cb_p_distrib = QtWidgets.QComboBox(self.groupBox)
        self.cb_p_distrib.setObjectName("cb_p_distrib")
        self.gridLayout.addWidget(self.cb_p_distrib, 0, 0, 1, 1)
        self.gridLayout_2.addWidget(self.groupBox, 1, 3, 1, 1)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout.setSpacing(6)
        self.verticalLayout.setObjectName("verticalLayout")
        self.cb_p_file = QtWidgets.QComboBox(Widget)
        self.cb_p_file.setObjectName("cb_p_file")
        self.verticalLayout.addWidget(self.cb_p_file)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout.setSpacing(6)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.pb_plot = QtWidgets.QPushButton(Widget)
        self.pb_plot.setStyleSheet("color: rgb(85, 255, 127);")
        self.pb_plot.setObjectName("pb_plot")
        self.horizontalLayout.addWidget(self.pb_plot)
        self.pb_reload_files = QtWidgets.QPushButton(Widget)
        self.pb_reload_files.setObjectName("pb_reload_files")
        self.horizontalLayout.addWidget(self.pb_reload_files)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.gridLayout_2.addLayout(self.verticalLayout, 1, 1, 1, 2)

        self.retranslateUi(Widget)
        QtCore.QMetaObject.connectSlotsByName(Widget)

    def retranslateUi(self, Widget):
        _translate = QtCore.QCoreApplication.translate
        Widget.setWindowTitle(_translate("Widget", "Widget"))
        self.groupBox_2.setTitle(_translate("Widget", "Slice Params"))
        self.groupBox.setTitle(_translate("Widget", "Particle Distribution"))
        self.pb_plot.setText(_translate("Widget", "Plot"))
        self.pb_reload_files.setText(_translate("Widget", "Reload Files"))

