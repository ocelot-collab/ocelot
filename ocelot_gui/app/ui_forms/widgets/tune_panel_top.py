# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'app/ui_forms/widgets/tune_panel_top.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form_TunePanel_Top(object):
    def setupUi(self, Form_TunePanel_Top):
        Form_TunePanel_Top.setObjectName("Form_TunePanel_Top")
        Form_TunePanel_Top.resize(285, 70)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form_TunePanel_Top.sizePolicy().hasHeightForWidth())
        Form_TunePanel_Top.setSizePolicy(sizePolicy)
        Form_TunePanel_Top.setMinimumSize(QtCore.QSize(285, 70))
        Form_TunePanel_Top.setMaximumSize(QtCore.QSize(285, 70))
        Form_TunePanel_Top.setAutoFillBackground(False)
        Form_TunePanel_Top.setStyleSheet(".QWidget {background-color: rgb(238, 238, 236); border-top: 1px solid rgb(255, 255, 255);}")
        self.xbutton = QtWidgets.QPushButton(Form_TunePanel_Top)
        self.xbutton.setGeometry(QtCore.QRect(260, 5, 20, 23))
        self.xbutton.setObjectName("xbutton")
        self.check = QtWidgets.QCheckBox(Form_TunePanel_Top)
        self.check.setGeometry(QtCore.QRect(10, 11, 121, 21))
        self.check.setAutoFillBackground(False)
        self.check.setStyleSheet("background-color: rgb(238, 238, 236)")
        self.check.setText("")
        self.check.setObjectName("check")
        self.label_f = QtWidgets.QLabel(Form_TunePanel_Top)
        self.label_f.setGeometry(QtCore.QRect(140, 13, 45, 16))
        self.label_f.setStyleSheet("background-color: rgb(238, 238, 236)")
        self.label_f.setObjectName("label_f")
        self.factor = QtWidgets.QLineEdit(Form_TunePanel_Top)
        self.factor.setGeometry(QtCore.QRect(190, 10, 50, 23))
        self.factor.setObjectName("factor")
        self.label_c = QtWidgets.QLabel(Form_TunePanel_Top)
        self.label_c.setGeometry(QtCore.QRect(5, 45, 45, 15))
        self.label_c.setStyleSheet("background-color: rgb(238, 238, 236)")
        self.label_c.setObjectName("label_c")
        self.label_o = QtWidgets.QLabel(Form_TunePanel_Top)
        self.label_o.setGeometry(QtCore.QRect(142, 45, 65, 15))
        self.label_o.setStyleSheet("background-color: rgb(238, 238, 236)")
        self.label_o.setObjectName("label_o")
        self.curval = QtWidgets.QDoubleSpinBox(Form_TunePanel_Top)
        self.curval.setGeometry(QtCore.QRect(46, 40, 90, 24))
        self.curval.setDecimals(6)
        self.curval.setMinimum(-1000000000000000.0)
        self.curval.setMaximum(1e+16)
        self.curval.setSingleStep(0.01)
        self.curval.setObjectName("curval")
        self.oldval = QtWidgets.QLineEdit(Form_TunePanel_Top)
        self.oldval.setGeometry(QtCore.QRect(205, 40, 75, 23))
        self.oldval.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.oldval.setReadOnly(True)
        self.oldval.setObjectName("oldval")

        self.retranslateUi(Form_TunePanel_Top)
        QtCore.QMetaObject.connectSlotsByName(Form_TunePanel_Top)

    def retranslateUi(self, Form_TunePanel_Top):
        _translate = QtCore.QCoreApplication.translate
        Form_TunePanel_Top.setWindowTitle(_translate("Form_TunePanel_Top", "Form"))
        self.xbutton.setText(_translate("Form_TunePanel_Top", "x"))
        self.label_f.setText(_translate("Form_TunePanel_Top", "factor:"))
        self.factor.setText(_translate("Form_TunePanel_Top", "0.01"))
        self.label_c.setText(_translate("Form_TunePanel_Top", "value:"))
        self.label_o.setText(_translate("Form_TunePanel_Top", "old value:"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form_TunePanel_Top = QtWidgets.QWidget()
    ui = Ui_Form_TunePanel_Top()
    ui.setupUi(Form_TunePanel_Top)
    Form_TunePanel_Top.show()
    sys.exit(app.exec_())

