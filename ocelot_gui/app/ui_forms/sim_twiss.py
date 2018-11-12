# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'app/ui_forms/sim_twiss.ui'
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
        self.label.setGeometry(QtCore.QRect(10, 0, 121, 16))
        self.label.setObjectName("label")
        self.label_3 = QtWidgets.QLabel(Form)
        self.label_3.setGeometry(QtCore.QRect(880, 0, 154, 16))
        self.label_3.setObjectName("label_3")
        self.btn1 = QtWidgets.QPushButton(Form)
        self.btn1.setGeometry(QtCore.QRect(1014, 710, 80, 23))
        self.btn1.setObjectName("btn1")
        self.btn2 = QtWidgets.QPushButton(Form)
        self.btn2.setGeometry(QtCore.QRect(1100, 710, 80, 23))
        self.btn2.setObjectName("btn2")
        self.edit_tws_step = QtWidgets.QLineEdit(Form)
        self.edit_tws_step.setGeometry(QtCore.QRect(880, 20, 80, 23))
        self.edit_tws_step.setAlignment(QtCore.Qt.AlignCenter)
        self.edit_tws_step.setObjectName("edit_tws_step")
        self.tws_plot_widget = QtWidgets.QWidget(Form)
        self.tws_plot_widget.setGeometry(QtCore.QRect(10, 20, 850, 650))
        self.tws_plot_widget.setObjectName("tws_plot_widget")
        self.label_4 = QtWidgets.QLabel(Form)
        self.label_4.setGeometry(QtCore.QRect(880, 60, 321, 16))
        self.label_4.setObjectName("label_4")
        self.edit_delta = QtWidgets.QLineEdit(Form)
        self.edit_delta.setGeometry(QtCore.QRect(990, 80, 80, 23))
        self.edit_delta.setAlignment(QtCore.Qt.AlignCenter)
        self.edit_delta.setObjectName("edit_delta")
        self.plus_button = QtWidgets.QPushButton(Form)
        self.plus_button.setGeometry(QtCore.QRect(950, 80, 30, 23))
        self.plus_button.setObjectName("plus_button")
        self.minus_button = QtWidgets.QPushButton(Form)
        self.minus_button.setGeometry(QtCore.QRect(1080, 80, 30, 23))
        self.minus_button.setObjectName("minus_button")
        self.scroll = QtWidgets.QScrollArea(Form)
        self.scroll.setGeometry(QtCore.QRect(880, 110, 300, 560))
        self.scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
        self.scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scroll.setWidgetResizable(True)
        self.scroll.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.scroll.setObjectName("scroll")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 284, 558))
        self.scrollAreaWidgetContents.setAutoFillBackground(False)
        self.scrollAreaWidgetContents.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.verticalLayoutWidget = QtWidgets.QWidget(self.scrollAreaWidgetContents)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(0, 0, 281, 80))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.widgetArea = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.widgetArea.setContentsMargins(0, 0, 0, 2)
        self.widgetArea.setSpacing(0)
        self.widgetArea.setObjectName("widgetArea")
        self.scroll.setWidget(self.scrollAreaWidgetContents)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.label.setText(_translate("Form", "Twiss functions"))
        self.label_3.setText(_translate("Form", "Twiss function step, m"))
        self.btn1.setText(_translate("Form", "Update"))
        self.btn2.setText(_translate("Form", "Reset"))
        self.edit_tws_step.setText(_translate("Form", "0.0"))
        self.label_4.setText(_translate("Form", "Tuning elements list"))
        self.edit_delta.setText(_translate("Form", "0.0"))
        self.plus_button.setText(_translate("Form", "+"))
        self.minus_button.setText(_translate("Form", "-"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())

