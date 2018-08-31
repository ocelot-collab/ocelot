
from PyQt5 import QtGui, QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget

from app.parser import *

class EditField(QWidget):

    def __init__(self, MainWindow, parent=None):
        super().__init__()
        self.mw = MainWindow
        self.init_form()


    def init_form(self):
        
        label1 = QtGui.QLabel("Elements list")
        self.edit_elements = QtGui.QTextEdit()

        label2 = QtGui.QLabel("Sequences list (the last sequence is used as active lattice)")
        self.edit_cells = QtGui.QTextEdit()

        label3 = QtGui.QLabel("Number of superperiods")
        self.edit_nsuperperiods = QtGui.QLineEdit()
        self.edit_nsuperperiods.setFixedWidth(50)

        btn1 = QtGui.QPushButton("Update")
        btn1.clicked.connect(self.action_update_parameters)
        
        btn2 = QtGui.QPushButton("Cancel")
        btn2.clicked.connect(self.action_cancel)

        labelcb = QtGui.QLabel("Periodic solution")
        self.edit_periodic_solution = QtGui.QCheckBox()

        label4 = QtGui.QLabel("Initial Beam parameters")
        self.edit_beam = QtGui.QTextEdit()

        label5 = QtGui.QLabel("Initial Twiss parameters")
        self.edit_twiss = QtGui.QTextEdit()

        layout = QtWidgets.QGridLayout()

        layout.addWidget(label1, 0, 0, 1, 3)
        layout.addWidget(self.edit_elements, 1, 0, 3, 3)
        layout.addWidget(label2, 4, 0, 1, 3)
        layout.addWidget(self.edit_cells, 5, 0, 4, 3)
        
        layout.addWidget(labelcb, 5, 3, Qt.AlignBottom)
        layout.addWidget(self.edit_periodic_solution, 6, 3, Qt.AlignTop)
        layout.addWidget(label3, 7, 3, Qt.AlignBottom)
        layout.addWidget(self.edit_nsuperperiods, 8, 3, Qt.AlignTop)

        layout.addWidget(label4, 0, 3)
        layout.addWidget(self.edit_beam, 1, 3)

        layout.addWidget(label5, 2, 3)
        layout.addWidget(self.edit_twiss, 3, 3)

        # fix for format mesh
        for i in range(4):
            l = QtGui.QLabel('')
            layout.addWidget(l, 9, i)
            layout.setColumnMinimumWidth(i, 300)
            layout.setColumnStretch(i, 1)

        box_buttons = QtGui.QHBoxLayout()
        box_buttons.addWidget(btn1)
        box_buttons.addWidget(btn2)
        
        layout.addLayout(box_buttons, 10, 3)

        self.setLayout(layout)
        self.load_parameters()


    def load_parameters(self):
        """Load all parameters to window"""

        parcer = LatticeSave()

        # elements
        lines = parcer.parsing_elements(self.mw.lattice.elements)
        lines = lines[0:-2]
        self.edit_elements.setText(lines)

        # sequences
        lines = parcer.parsing_cells(self.mw.lattice.cells, self.mw.lattice.cells_order)
        lines = lines[0:-2]
        lines = lines.replace('\n', '\n\n')
        self.edit_cells.setText(lines)

        # number of superperiods
        self.edit_nsuperperiods.setText(str(self.mw.lattice.nsuperperiods))

        # beam
        lines = parcer.parsing_beam(self.mw.lattice.beam)
        lines = lines[0:-2]
        self.edit_beam.setText(lines)

        # tws
        lines = parcer.parsing_twiss(self.mw.lattice.tws)
        lines = lines[0:-2]
        self.edit_twiss.setText(lines)
        
        # periodic solution
        if self.mw.lattice.periodic_solution:
            self.edit_periodic_solution.toggle()


    def action_update_parameters(self):
        
        parcer = LatticeRead()

        # elements
        lines = parcer.del_comments(self.edit_elements.toPlainText())
        self.mw.lattice.elements = parcer.parsing_elements(lines)

        # sequences
        lines = parcer.del_comments(self.edit_cells.toPlainText())
        self.mw.lattice.cells, self.mw.lattice.cells_order = parcer.parsing_cells(lines, self.mw.lattice.elements)

        # beam
        lines = parcer.del_comments(self.edit_beam.toPlainText())
        self.mw.lattice.beam = parcer.parsing_beam(lines)

        # tws
        lines = parcer.del_comments(self.edit_twiss.toPlainText())
        self.mw.lattice.tws = parcer.parsing_twiss(lines)

        # number of superperiods
        try:
            val = int(self.edit_nsuperperiods.text())
            if val > 0:
                self.mw.lattice.nsuperperiods = val
        except ValueError:
            self.mw.lattice.nsuperperiods = 1

        # periodic solution
        if self.edit_periodic_solution.checkState():
            self.mw.lattice.periodic_solution = True
        else:
            self.mw.lattice.periodic_solution = False

        self.mw.init_central_widget()

    
    def action_cancel(self):
        self.mw.init_central_widget()
