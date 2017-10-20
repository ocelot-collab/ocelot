# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'UIresetpanel.ui'
#
# Created: Wed Jun 15 13:36:17 2016
#      by: PyQt4 UI code generator 4.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(500, 850)
        Form.setMinimumSize(QtCore.QSize(400, 0))
        Form.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        Form.setStyleSheet(_fromUtf8("background-color: white"))
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.tableWidget = QtGui.QTableWidget(Form)
        self.tableWidget.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.tableWidget.setObjectName(_fromUtf8("tableWidget"))
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)
        self.gridLayout.addWidget(self.tableWidget, 0, 0, 1, 1)
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.updateReference = QtGui.QPushButton(Form)
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("DejaVu Sans"))
        font.setPointSize(10)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.updateReference.setFont(font)
        self.updateReference.setStyleSheet(_fromUtf8("color: orange"))
        self.updateReference.setObjectName(_fromUtf8("updateReference"))
        self.horizontalLayout.addWidget(self.updateReference)
        self.resetAll = QtGui.QPushButton(Form)
        font = QtGui.QFont()
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.resetAll.setFont(font)
        self.resetAll.setStyleSheet(_fromUtf8("color: red"))
        self.resetAll.setObjectName(_fromUtf8("resetAll"))
        self.horizontalLayout.addWidget(self.resetAll)
        self.gridLayout.addLayout(self.horizontalLayout, 2, 0, 1, 1)
        self.label = QtGui.QLabel(Form)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setItalic(True)
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 1, 0, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Reset Panel", None))
        self.updateReference.setText(_translate("Form", "Update reference ", None))
        self.resetAll.setText(_translate("Form", "Reset All", None))
        self.label.setText(_translate("Form", "Middle click a PV then the table to add your favorite device!", None))

