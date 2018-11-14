# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'tune_panel_bottom.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form_TunePanel_Bottom(object):
    def setupUi(self, Form_TunePanel_Bottom):
        Form_TunePanel_Bottom.setObjectName("Form_TunePanel_Bottom")
        Form_TunePanel_Bottom.resize(285, 70)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form_TunePanel_Bottom.sizePolicy().hasHeightForWidth())
        Form_TunePanel_Bottom.setSizePolicy(sizePolicy)
        Form_TunePanel_Bottom.setMinimumSize(QtCore.QSize(285, 70))
        Form_TunePanel_Bottom.setMaximumSize(QtCore.QSize(285, 70))
        Form_TunePanel_Bottom.setAutoFillBackground(False)
        Form_TunePanel_Bottom.setStyleSheet(".QWidget {background-color: rgb(238, 238, 236);}")
        self.gridLayout = QtWidgets.QGridLayout(Form_TunePanel_Bottom)
        self.gridLayout.setObjectName("gridLayout")
        self.check = QtWidgets.QCheckBox(Form_TunePanel_Bottom)
        self.check.setStyleSheet("background-color: rgb(238, 238, 236)")
        self.check.setText("")
        self.check.setObjectName("check")
        self.gridLayout.addWidget(self.check, 0, 0, 1, 2)
        self.label_f = QtWidgets.QLabel(Form_TunePanel_Bottom)
        self.label_f.setStyleSheet("background-color: rgb(238, 238, 236)")
        self.label_f.setObjectName("label_f")
        self.gridLayout.addWidget(self.label_f, 0, 2, 1, 1)
        self.factor = QtWidgets.QLineEdit(Form_TunePanel_Bottom)
        self.factor.setMinimumSize(QtCore.QSize(50, 23))
        self.factor.setObjectName("factor")
        self.gridLayout.addWidget(self.factor, 0, 3, 1, 2)
        self.label_c = QtWidgets.QLabel(Form_TunePanel_Bottom)
        self.label_c.setStyleSheet("background-color: rgb(238, 238, 236)")
        self.label_c.setObjectName("label_c")
        self.gridLayout.addWidget(self.label_c, 1, 0, 1, 1)
        self.curval = QtWidgets.QDoubleSpinBox(Form_TunePanel_Bottom)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.curval.sizePolicy().hasHeightForWidth())
        self.curval.setSizePolicy(sizePolicy)
        self.curval.setMinimumSize(QtCore.QSize(0, 23))
        self.curval.setMaximumSize(QtCore.QSize(95, 16777215))
        self.curval.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
        self.curval.setDecimals(6)
        self.curval.setMinimum(-1000000000000000.0)
        self.curval.setMaximum(1e+16)
        self.curval.setSingleStep(0.01)
        self.curval.setObjectName("curval")
        self.gridLayout.addWidget(self.curval, 1, 1, 1, 1)
        self.label_o = QtWidgets.QLabel(Form_TunePanel_Bottom)
        self.label_o.setStyleSheet("background-color: rgb(238, 238, 236)")
        self.label_o.setObjectName("label_o")
        self.gridLayout.addWidget(self.label_o, 1, 2, 1, 2)
        self.oldval = QtWidgets.QLineEdit(Form_TunePanel_Bottom)
        self.oldval.setMinimumSize(QtCore.QSize(0, 23))
        self.oldval.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.oldval.setReadOnly(True)
        self.oldval.setObjectName("oldval")
        self.gridLayout.addWidget(self.oldval, 1, 4, 1, 1)

        self.retranslateUi(Form_TunePanel_Bottom)
        QtCore.QMetaObject.connectSlotsByName(Form_TunePanel_Bottom)

    def retranslateUi(self, Form_TunePanel_Bottom):
        _translate = QtCore.QCoreApplication.translate
        Form_TunePanel_Bottom.setWindowTitle(_translate("Form_TunePanel_Bottom", "Form"))
        self.label_f.setText(_translate("Form_TunePanel_Bottom", "factor:"))
        self.factor.setText(_translate("Form_TunePanel_Bottom", "0.01"))
        self.label_c.setText(_translate("Form_TunePanel_Bottom", "value:"))
        self.label_o.setText(_translate("Form_TunePanel_Bottom", "old value:"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form_TunePanel_Bottom = QtWidgets.QWidget()
    ui = Ui_Form_TunePanel_Bottom()
    ui.setupUi(Form_TunePanel_Bottom)
    Form_TunePanel_Bottom.show()
    sys.exit(app.exec_())

