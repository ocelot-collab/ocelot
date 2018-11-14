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
        Form_MatchingPanel.resize(285, 35)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form_MatchingPanel.sizePolicy().hasHeightForWidth())
        Form_MatchingPanel.setSizePolicy(sizePolicy)
        Form_MatchingPanel.setMinimumSize(QtCore.QSize(285, 35))
        Form_MatchingPanel.setMaximumSize(QtCore.QSize(285, 35))
        Form_MatchingPanel.setAutoFillBackground(False)
        Form_MatchingPanel.setStyleSheet(".QWidget {background-color: rgb(238, 238, 236); border-top: 1px solid rgb(255, 255, 255);}")
        self.horizontalLayoutWidget = QtWidgets.QWidget(Form_MatchingPanel)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(0, 0, 282, 31))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout_2.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
        self.horizontalLayout_2.setContentsMargins(0, 5, 0, 0)
        self.horizontalLayout_2.setSpacing(0)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        spacerItem = QtWidgets.QSpacerItem(5, 23, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.check = QtWidgets.QCheckBox(self.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.check.sizePolicy().hasHeightForWidth())
        self.check.setSizePolicy(sizePolicy)
        self.check.setMinimumSize(QtCore.QSize(120, 23))
        self.check.setMaximumSize(QtCore.QSize(120, 23))
        self.check.setAutoFillBackground(False)
        self.check.setStyleSheet("background-color: rgb(238, 238, 236)")
        self.check.setText("")
        self.check.setChecked(True)
        self.check.setObjectName("check")
        self.horizontalLayout_2.addWidget(self.check)
        spacerItem1 = QtWidgets.QSpacerItem(8, 23, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.label_c = QtWidgets.QLabel(self.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_c.sizePolicy().hasHeightForWidth())
        self.label_c.setSizePolicy(sizePolicy)
        self.label_c.setMinimumSize(QtCore.QSize(40, 23))
        self.label_c.setMaximumSize(QtCore.QSize(40, 23))
        self.label_c.setStyleSheet("background-color: rgb(238, 238, 236)")
        self.label_c.setObjectName("label_c")
        self.horizontalLayout_2.addWidget(self.label_c)
        self.oldval = QtWidgets.QLineEdit(self.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.oldval.sizePolicy().hasHeightForWidth())
        self.oldval.setSizePolicy(sizePolicy)
        self.oldval.setMinimumSize(QtCore.QSize(75, 23))
        self.oldval.setMaximumSize(QtCore.QSize(75, 23))
        self.oldval.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.oldval.setReadOnly(True)
        self.oldval.setObjectName("oldval")
        self.horizontalLayout_2.addWidget(self.oldval)
        spacerItem2 = QtWidgets.QSpacerItem(9, 23, QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem2)
        self.xbutton = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.xbutton.sizePolicy().hasHeightForWidth())
        self.xbutton.setSizePolicy(sizePolicy)
        self.xbutton.setMinimumSize(QtCore.QSize(23, 23))
        self.xbutton.setMaximumSize(QtCore.QSize(23, 23))
        self.xbutton.setObjectName("xbutton")
        self.horizontalLayout_2.addWidget(self.xbutton)
        self.horizontalLayoutWidget.raise_()
        self.check.raise_()

        self.retranslateUi(Form_MatchingPanel)
        QtCore.QMetaObject.connectSlotsByName(Form_MatchingPanel)

    def retranslateUi(self, Form_MatchingPanel):
        _translate = QtCore.QCoreApplication.translate
        Form_MatchingPanel.setWindowTitle(_translate("Form_MatchingPanel", "Form"))
        self.label_c.setText(_translate("Form_MatchingPanel", "value:"))
        self.xbutton.setText(_translate("Form_MatchingPanel", "x"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form_MatchingPanel = QtWidgets.QWidget()
    ui = Ui_Form_MatchingPanel()
    ui.setupUi(Form_MatchingPanel)
    Form_MatchingPanel.show()
    sys.exit(app.exec_())

