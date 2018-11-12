# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'app/forms/ui_forms/editwindow.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(1200, 750)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        self.label = QtWidgets.QLabel(Form)
        self.label.setGeometry(QtCore.QRect(10, 0, 82, 16))
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(Form)
        self.label_2.setGeometry(QtCore.QRect(10, 515, 86, 16))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(Form)
        self.label_3.setGeometry(QtCore.QRect(880, 0, 154, 16))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(Form)
        self.label_4.setGeometry(QtCore.QRect(880, 240, 149, 16))
        self.label_4.setObjectName("label_4")
        self.edit_periodic_solution = QtWidgets.QCheckBox(Form)
        self.edit_periodic_solution.setGeometry(QtCore.QRect(881, 550, 129, 21))
        self.edit_periodic_solution.setObjectName("edit_periodic_solution")
        self.label_5 = QtWidgets.QLabel(Form)
        self.label_5.setGeometry(QtCore.QRect(881, 610, 155, 16))
        self.label_5.setObjectName("label_5")
        self.edit_nsuperperiods = QtWidgets.QSpinBox(Form)
        self.edit_nsuperperiods.setGeometry(QtCore.QRect(881, 630, 80, 24))
        self.edit_nsuperperiods.setMaximumSize(QtCore.QSize(80, 16777215))
        self.edit_nsuperperiods.setMinimum(1)
        self.edit_nsuperperiods.setMaximum(999999)
        self.edit_nsuperperiods.setObjectName("edit_nsuperperiods")
        self.btn1 = QtWidgets.QPushButton(Form)
        self.btn1.setGeometry(QtCore.QRect(1014, 710, 80, 23))
        self.btn1.setObjectName("btn1")
        self.edit_elements = QtWidgets.QTextEdit(Form)
        self.edit_elements.setGeometry(QtCore.QRect(10, 20, 850, 491))
        self.edit_elements.setObjectName("edit_elements")
        self.edit_cells = QtWidgets.QTextEdit(Form)
        self.edit_cells.setGeometry(QtCore.QRect(10, 535, 850, 192))
        self.edit_cells.setObjectName("edit_cells")
        self.edit_beam = QtWidgets.QTextEdit(Form)
        self.edit_beam.setGeometry(QtCore.QRect(880, 20, 300, 210))
        self.edit_beam.setObjectName("edit_beam")
        self.edit_twiss = QtWidgets.QTextEdit(Form)
        self.edit_twiss.setGeometry(QtCore.QRect(880, 260, 300, 250))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.edit_twiss.sizePolicy().hasHeightForWidth())
        self.edit_twiss.setSizePolicy(sizePolicy)
        self.edit_twiss.setObjectName("edit_twiss")
        self.btn2 = QtWidgets.QPushButton(Form)
        self.btn2.setGeometry(QtCore.QRect(1100, 710, 80, 23))
        self.btn2.setObjectName("btn2")
        self.label.raise_()
        self.label_2.raise_()
        self.label_3.raise_()
        self.label_4.raise_()
        self.edit_periodic_solution.raise_()
        self.label_5.raise_()
        self.edit_nsuperperiods.raise_()
        self.btn1.raise_()
        self.btn2.raise_()
        self.edit_elements.raise_()
        self.edit_cells.raise_()
        self.edit_beam.raise_()
        self.edit_twiss.raise_()
        self.btn2.raise_()

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.label.setText(_translate("Form", "Elements list"))
        self.label_2.setText(_translate("Form", "Sequence list"))
        self.label_3.setText(_translate("Form", "Initial Beam parameters"))
        self.label_4.setText(_translate("Form", "Initial Twiss parameters"))
        self.edit_periodic_solution.setText(_translate("Form", "Periodic solution"))
        self.label_5.setText(_translate("Form", "Number of superperiods"))
        self.btn1.setText(_translate("Form", "Update"))
        self.btn2.setText(_translate("Form", "Cancel"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())

