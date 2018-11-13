# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'app/ui_forms/widgets/matching_panel.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form_MatchingPanel(object):
    def setupUi(self, Form_MatchingPanel):
        Form_MatchingPanel.setObjectName("Form_MatchingPanel")
        Form_MatchingPanel.resize(285, 45)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form_MatchingPanel.sizePolicy().hasHeightForWidth())
        Form_MatchingPanel.setSizePolicy(sizePolicy)
        Form_MatchingPanel.setMinimumSize(QtCore.QSize(285, 45))
        Form_MatchingPanel.setMaximumSize(QtCore.QSize(285, 45))
        Form_MatchingPanel.setAutoFillBackground(False)
        Form_MatchingPanel.setStyleSheet(".QWidget {background-color: rgb(238, 238, 236); border-top: 1px solid rgb(255, 255, 255);}")
        self.xbutton = QtWidgets.QPushButton(Form_MatchingPanel)
        self.xbutton.setGeometry(QtCore.QRect(260, 5, 20, 23))
        self.xbutton.setObjectName("xbutton")
        self.check = QtWidgets.QCheckBox(Form_MatchingPanel)
        self.check.setGeometry(QtCore.QRect(10, 11, 121, 21))
        self.check.setAutoFillBackground(False)
        self.check.setStyleSheet("background-color: rgb(238, 238, 236)")
        self.check.setText("")
        self.check.setChecked(True)
        self.check.setObjectName("check")
        self.label_c = QtWidgets.QLabel(Form_MatchingPanel)
        self.label_c.setGeometry(QtCore.QRect(140, 15, 45, 15))
        self.label_c.setStyleSheet("background-color: rgb(238, 238, 236)")
        self.label_c.setObjectName("label_c")
        self.oldval = QtWidgets.QLineEdit(Form_MatchingPanel)
        self.oldval.setGeometry(QtCore.QRect(180, 12, 75, 23))
        self.oldval.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.oldval.setReadOnly(True)
        self.oldval.setObjectName("oldval")

        self.retranslateUi(Form_MatchingPanel)
        QtCore.QMetaObject.connectSlotsByName(Form_MatchingPanel)

    def retranslateUi(self, Form_MatchingPanel):
        _translate = QtCore.QCoreApplication.translate
        Form_MatchingPanel.setWindowTitle(_translate("Form_MatchingPanel", "Form"))
        self.xbutton.setText(_translate("Form_MatchingPanel", "x"))
        self.label_c.setText(_translate("Form_MatchingPanel", "value:"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form_MatchingPanel = QtWidgets.QWidget()
    ui = Ui_Form_MatchingPanel()
    ui.setupUi(Form_MatchingPanel)
    Form_MatchingPanel.show()
    sys.exit(app.exec_())

