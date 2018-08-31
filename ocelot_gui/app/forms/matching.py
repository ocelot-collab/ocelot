from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QWidget, QTextEdit, QLabel, QPushButton

from copy import deepcopy

from ocelot import *

from app.forms.widgets.lattice_plot import *
from app.forms.widgets.vars_panel import *
from app.forms.widgets.constr_panel import *


class Matching(QWidget):

    def __init__(self, MainWindow, parent=None):
        super().__init__()
        self.mw = MainWindow
        self.vars_elements = {}
        self.init_form()


    def init_form(self):

        self.mw.lattice.init_lattice()

        self.vars_panel = VarsPanel(self.mw, self.vars_elements)
        self.lattice_plot = LatticePlot(self.mw, self.vars_panel)
        self.constraints = ConstraintsPanel(self.mw)
        
        label6 = QLabel("Matching results")
        self.results = QTextEdit()
        self.results.setMinimumHeight(250)
        self.results.setReadOnly(True)

        btn1 = QPushButton("Apply results")
        btn1.clicked.connect(self.action_apply_results)

        btn2 = QPushButton("Matching")
        btn2.clicked.connect(self.action_match)

        layout = QtGui.QGridLayout()
        layout.addWidget(self.lattice_plot.lattice_plot, 0, 0, 1, 3)
        layout.addLayout(self.vars_panel.vars_panel, 0, 3, 4, 1)
        layout.addLayout(self.constraints.constr_panel, 2, 0, 1, 3)

        layout.addWidget(label6, 3, 0)
        layout.addWidget(self.results, 4, 0, 2, 3)

        # fix for format mesh
        for i in range(4):
            l = QtGui.QLabel('')
            layout.addWidget(l, 6, i)
            layout.setColumnMinimumWidth(i, 300)
            layout.setColumnStretch(i, 1)

        box_buttons = QtGui.QHBoxLayout()
        box_buttons.addWidget(btn1)
        box_buttons.addWidget(btn2)
        
        layout.addLayout(box_buttons, 7, 3)

        self.setLayout(layout)


    def action_match(self):

        self.mw.work_lattice = deepcopy(self.mw.lattice)

        tws = Twiss()
        tws.E = self.mw.lattice.tws.E
        tws_labels = ['beta_x', 'beta_y', 'alpha_x', 'alpha_y', 'Dx', 'Dxp']
        m_start = Monitor(eid="matching_start")
        m_end = Monitor(eid="matching_end")

        # collect variables
        vars = []
        for key in self.vars_elements:
            name = self.vars_elements[key]['element']
            vars.append(self.mw.work_lattice.elements[name])

        if vars == []:
            self.results.setText('Select matching variables')
            return

        # collect constraints
        constr = {}

        # check periodic constraint
        if self.constraints.periodic_solution.checkState():
            constr['periodic'] = True

        # check global constraints
        g = {}
        for i in range(6):
            val = self.constraints.twiss[2][i].text()
            if val != '':
                g[tws_labels[i]] = [str(self.constraints.conds[2][i].currentText()), float(val)]
        
        if g!= {}:
            constr['global'] = g

        # check constraints at lattice start and end
        g1, g2 = {}, {}
        for i in range(6):

            val1 = self.constraints.twiss[0][i].text()
            val2 = self.constraints.twiss[1][i].text()

            cond1 = str(self.constraints.conds[0][i].currentText())
            cond2 = str(self.constraints.conds[1][i].currentText())
            
            if val1 != '':
                tws.__dict__[tws_labels[i]] = float(val1)
                g1[tws_labels[i]] = float(val1) if cond1 == '=' else [str(self.constraints.conds[0][i].currentText()), float(val1)]
            
            if val2 != '':
                g2[tws_labels[i]] = float(val2) if cond2 == '=' else [str(self.constraints.conds[1][i].currentText()), float(val2)]

        if g1!= {}:
            constr[m_start] = g1
        if g2!= {}:
            constr[m_end] = g2

        if constr == {}:
            self.results.setText('Set matching constraints')
            return

        # prepare lattice
        self.mw.work_lattice.tws = tws
        self.mw.work_lattice.init_lattice()

        self.mw.work_lattice.elements[m_start.id] = m_start
        self.mw.work_lattice.elements[m_end.id] = m_end

        sec_name = self.mw.work_lattice.cells_order[-1]
        self.mw.work_lattice.cells[sec_name] = (m_start.id,) + self.mw.work_lattice.cells[sec_name] + (m_end.id,)
        self.mw.work_lattice.init_lattice()

        # matching
        result = match(self.mw.work_lattice.lattice, constr, vars, tws)

        lines = 'New variable values:\n'
        for i, v in enumerate(vars):

            eclass = v.__class__.__name__
            type = self.mw.matchable_elements[eclass][0]
            lines += str(v.id) + '.' + str(type) + ' = ' + str(result[i]) + '\n'

        self.results.setText(lines)


    def action_apply_results(self):
        
        self.mw.lattice = deepcopy(self.mw.work_lattice)
        self.results.setText('Lattice has updated')
