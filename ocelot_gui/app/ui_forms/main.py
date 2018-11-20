# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'ocelot_gui/app/ui_forms/main.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1200, 800)
        MainWindow.setMinimumSize(QtCore.QSize(1150, 0))
        self.central_widget = QtWidgets.QWidget(MainWindow)
        self.central_widget.setObjectName("central_widget")
        MainWindow.setCentralWidget(self.central_widget)
        self.statusBar = QtWidgets.QStatusBar(MainWindow)
        self.statusBar.setObjectName("statusBar")
        MainWindow.setStatusBar(self.statusBar)
        self.menuBar = QtWidgets.QMenuBar(MainWindow)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 1200, 20))
        self.menuBar.setObjectName("menuBar")
        self.menuFile = QtWidgets.QMenu(self.menuBar)
        self.menuFile.setObjectName("menuFile")
        self.menuEdit = QtWidgets.QMenu(self.menuBar)
        self.menuEdit.setObjectName("menuEdit")
        self.menuSimulations = QtWidgets.QMenu(self.menuBar)
        self.menuSimulations.setObjectName("menuSimulations")
        MainWindow.setMenuBar(self.menuBar)
        self.action_new_lattice = QtWidgets.QAction(MainWindow)
        self.action_new_lattice.setObjectName("action_new_lattice")
        self.action_open_lattice = QtWidgets.QAction(MainWindow)
        self.action_open_lattice.setObjectName("action_open_lattice")
        self.action_save_lattice = QtWidgets.QAction(MainWindow)
        self.action_save_lattice.setObjectName("action_save_lattice")
        self.action_exit = QtWidgets.QAction(MainWindow)
        self.action_exit.setObjectName("action_exit")
        self.action_edit_lattice = QtWidgets.QAction(MainWindow)
        self.action_edit_lattice.setObjectName("action_edit_lattice")
        self.action_calc_twiss = QtWidgets.QAction(MainWindow)
        self.action_calc_twiss.setObjectName("action_calc_twiss")
        self.action_calc_matching = QtWidgets.QAction(MainWindow)
        self.action_calc_matching.setObjectName("action_calc_matching")
        self.action_calc_params = QtWidgets.QAction(MainWindow)
        self.action_calc_params.setObjectName("action_calc_params")
        self.menuFile.addAction(self.action_new_lattice)
        self.menuFile.addAction(self.action_open_lattice)
        self.menuFile.addAction(self.action_save_lattice)
        self.menuFile.addAction(self.action_exit)
        self.menuEdit.addAction(self.action_edit_lattice)
        self.menuSimulations.addAction(self.action_calc_twiss)
        self.menuSimulations.addAction(self.action_calc_params)
        self.menuSimulations.addAction(self.action_calc_matching)
        self.menuBar.addAction(self.menuFile.menuAction())
        self.menuBar.addAction(self.menuEdit.menuAction())
        self.menuBar.addAction(self.menuSimulations.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuEdit.setTitle(_translate("MainWindow", "Edit"))
        self.menuSimulations.setTitle(_translate("MainWindow", "Simulations"))
        self.action_new_lattice.setText(_translate("MainWindow", "New"))
        self.action_open_lattice.setText(_translate("MainWindow", "Open Lattice"))
        self.action_save_lattice.setText(_translate("MainWindow", "Save Lattice"))
        self.action_exit.setText(_translate("MainWindow", "Exit"))
        self.action_edit_lattice.setText(_translate("MainWindow", "Edit Lattice"))
        self.action_edit_lattice.setToolTip(_translate("MainWindow", "Edit Lattice and parameters"))
        self.action_calc_twiss.setText(_translate("MainWindow", "Twiss Functions"))
        self.action_calc_twiss.setToolTip(_translate("MainWindow", "Twiss Functions"))
        self.action_calc_matching.setText(_translate("MainWindow", "Matching"))
        self.action_calc_params.setText(_translate("MainWindow", "Main Parameters"))
        self.action_calc_params.setToolTip(_translate("MainWindow", "Main Parameters Calculation"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

