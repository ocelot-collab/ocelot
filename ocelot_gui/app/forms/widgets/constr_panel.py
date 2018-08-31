from PyQt5 import QtGui
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget
from functools import partial

from copy import deepcopy


class ConstraintsPanel():

    def __init__(self, MainWindow):
        self.mw = MainWindow
        self.init_panel()


    def init_panel(self):
        
        self.constr_panel = QtGui.QVBoxLayout()
        self.constr_panel.setAlignment(Qt.AlignTop)

        label_title = QtGui.QLabel('Constraints')
        self.constr_panel.addWidget(label_title)

        self.layout = QtGui.QGridLayout()

        labelcb = QtGui.QLabel("Periodic solution")
        labelcb.setFixedWidth(150)
        self.periodic_solution = QtGui.QCheckBox()

        titles = []
        titles.append([QtGui.QLabel("At lattice start"), ["=", "<", ">"]])
        titles.append([QtGui.QLabel("At lattice end"),["=", "<", ">"]])
        titles.append([QtGui.QLabel("Along all lattice"), ["<", ">"]])

        self.conds, self.twiss = [], []
        
        for j in range(3):
           
            c_arr, t_arr = [], []
            titles[j][0].setFixedWidth(150)

            for i in range(6):
                
                c_arr.append(QtGui.QComboBox())
                c_arr[i].addItems(titles[j][1])
                c_arr[i].setFixedWidth(45)

                t_arr.append(QtGui.QLineEdit())
                t_arr[i].setFixedWidth(75)

            self.twiss.append(t_arr)
            self.conds.append(c_arr)

        self.layout.addWidget(labelcb, 0, 0)
        self.layout.addWidget(self.periodic_solution, 1, 0)
        
        for i in range(3):
            self.prepare_columbs(titles[i][0], i+1)

        self.constr_panel.addLayout(self.layout)


    def prepare_columbs(self, title, pos):
            
        labels_f = []
        labels_f.append(QtGui.QLabel("beta_x"))
        labels_f.append(QtGui.QLabel("beta_y"))
        labels_f.append(QtGui.QLabel("alpha_x"))
        labels_f.append(QtGui.QLabel("alpha_y"))
        labels_f.append(QtGui.QLabel("Dx"))
        labels_f.append(QtGui.QLabel("Dxp"))

        for i in range(6):
            labels_f[i].setFixedWidth(50)
            self.layout.addWidget(title, 0, pos)

            if pos == 1:
                val = self.mw.lattice.tws.__dict__[labels_f[i].text()]
                if val != 0.0:
                   self.twiss[pos-1][i].setText(str(round(val, 6)))

            hbox = QtGui.QHBoxLayout()
            hbox.addStretch(1)
            hbox.addWidget(labels_f[i])
            hbox.addWidget(self.conds[pos-1][i])
            hbox.addWidget(self.twiss[pos-1][i])
            hbox.addStretch(3)
            
            self.layout.addLayout(hbox, i+1, pos)
