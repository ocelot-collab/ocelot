from PyQt5 import QtGui
from PyQt5.QtWidgets import QWidget

from copy import deepcopy

from app.forms.widgets.twiss_plot import *
from app.forms.widgets.ebeam_table import *
from app.forms.widgets.tune_panel import *


class TwissMonitor(QWidget):

    def __init__(self, MainWindow, parent=None):
        super().__init__()
        self.mw = MainWindow
        self.tune_elements = {}
        self.tws_split = 0.1
        self.mw.work_lattice = deepcopy(self.mw.lattice)
        self.funcs = [self.change_factor, self.change_value, self.reset_lattice, self.update_lattice, self.change_value_up_down]
        self.init_form()


    def init_form(self):
        
        self.mw.work_lattice.init_lattice()

        self.calc_twiss()
        self.calc_ebeam_params()

        self.table = EbeamTable()
        self.tune_panel = TunePanel(self.mw, self.tune_elements, self.funcs)
        self.tws_plot = TwissPlot(self.mw, self.tune_panel)

        self.tws_plot.update_plot(self.tws)
        self.table.update_table(self.ebp)

        layout = QtGui.QGridLayout()
        layout.addWidget(self.tws_plot.tws_plot, 0, 0, 2, 3)
        layout.addWidget(self.table.table, 0, 3)
        layout.addLayout(self.tune_panel.tune_panel, 1, 3)

        # fix for format mesh
        for i in range(4):
            l = QtGui.QLabel('')
            layout.addWidget(l, 2, i)
            layout.setColumnMinimumWidth(i, 300)
            layout.setColumnStretch(i, 1)

        self.setLayout(layout)


    def calc_twiss(self):

        if self.mw.work_lattice.periodic_solution:
            tws0 = periodic_twiss(Twiss(), lattice_transfer_map(self.mw.work_lattice.lattice, self.mw.work_lattice.beam.E))
        else:
            tws0 = self.mw.work_lattice.tws

        if tws0 is not None:
            n = int(self.mw.work_lattice.lattice.totalLen / self.tws_split)
            self.tws = twiss(self.mw.work_lattice.lattice, tws0=tws0, nPoints=n)
        else:
            self.tws = None


    def init_ebp(self):

        self.ebp = []
        self.ebp.append([None])
        self.ebp.append(['Energy, GeV', '', 2])
        self.ebp.append(['Lorentz factor', '', 2])
        self.ebp.append(['Periods', '', 0])
        self.ebp.append(['Tolal length, m', '', 3])
        self.ebp.append(['Emittance, pm*rad', '', 2])
        self.ebp.append(['Tune X / Y', '', 4])
        self.ebp.append(['Natural chromaticity', '', 2])
        self.ebp.append(['Alpha, 1.0e-3', '', 3])
        self.ebp.append(['Energy spread, 1.0e-3', '', 3])
        self.ebp.append(['Energy loss, MeV', '', 3])
        self.ebp.append(['Jx / Jy / Je', '', 3])
        self.ebp.append(['tau X / Y / E ms', '', 3])


    def calc_ebeam_params(self):

        self.init_ebp()

        if self.tws is None:
            self.ebp[0] = None
            return
        
        (I1, I2, I3, I4, I5) = radiation_integrals(self.mw.work_lattice.lattice, self.tws[0], self.mw.work_lattice.nsuperperiods)
        if I2 == 0.0:
            self.ebp[0] = None
            return

        Q = [self.tws[-1].mux/2.0/np.pi*self.mw.work_lattice.nsuperperiods, self.tws[-1].muy/2.0/np.pi*self.mw.work_lattice.nsuperperiods]
        q_ksi = natural_chromaticity(self.mw.work_lattice.lattice, self.tws[0], self.mw.work_lattice.nsuperperiods)
        alpha = I1 / (self.mw.work_lattice.nsuperperiods * self.mw.work_lattice.lattice.totalLen) * 1.0e3
        Jx, Jy, Je = 1.0 - I4 / I2, 1.0, 2.0 + I4 / I2

        self.ebp[3][1] = self.mw.work_lattice.nsuperperiods
        self.ebp[4][1] = round(self.mw.work_lattice.nsuperperiods*self.mw.work_lattice.lattice.totalLen, self.ebp[4][2])
        self.ebp[6][1] = str(round(Q[0], self.ebp[6][2])) + ' / ' + str(round(Q[1], self.ebp[6][2]))
        self.ebp[7][1] = str(round(q_ksi[0], self.ebp[7][2])) + ' / ' + str(round(q_ksi[1], self.ebp[7][2]))
        self.ebp[8][1] = round(alpha, self.ebp[8][2])
        self.ebp[11][1] = str(round(Jx, self.ebp[11][2])) + ' / ' + str(round(Jy, self.ebp[11][2])) + ' / ' + str(round(Je, self.ebp[11][2]))

        if self.mw.work_lattice.beam.E > 0.0:
            
            gamma = self.mw.work_lattice.beam.E/m_e_GeV
            emittance = Cq * gamma**2 * I5 / Jx / I2 * 1.0e12
            u0 = Cgamma * (self.mw.work_lattice.beam.E * 1000.0)**4 * I2 / 2.0 / np.pi
            sigma_e = gamma * np.sqrt(Cq * I3 / Je / I2) * 1.0e3
            tau0 = 2.0 * self.mw.work_lattice.beam.E * 1000.0 * self.mw.work_lattice.nsuperperiods * self.mw.work_lattice.lattice.totalLen / speed_of_light / u0 * 1.0e3

            self.ebp[1][1] = round(self.mw.work_lattice.beam.E, self.ebp[1][2])
            self.ebp[2][1] = round(gamma, self.ebp[12][2])
            self.ebp[5][1] = round(emittance, self.ebp[5][2])
            self.ebp[9][1] = round(sigma_e, self.ebp[9][2])
            self.ebp[10][1] = round(u0, self.ebp[10][2])
            self.ebp[12][1] = str(round(tau0/Jx, self.ebp[12][2])) + ' / ' + str(round(tau0/Jy, self.ebp[12][2])) + ' / ' + str(round(tau0/Je, self.ebp[12][2]))    


    def change_factor(self, id):
        
        step = float(self.tune_elements[id[0]][id[1]]['factor'].text())
        self.tune_elements[id[0]][id[1]]['val'].setSingleStep(step)

    
    def change_value(self, id, d):

        self.tune_elements[id[0]]['element'].__dict__[id[1]] = float(d)
        self.mw.work_lattice.lattice.update_transfer_maps()
        
        self.calc_twiss()
        self.tws_plot.update_plot(self.tws)

        self.calc_ebeam_params()
        self.table.update_table(self.ebp)


    def reset_lattice(self):

        self.mw.work_lattice = deepcopy(self.mw.lattice)
        self.mw.work_lattice.init_lattice()

        self.calc_twiss()
        self.tws_plot.update_plot(self.tws)

        self.calc_ebeam_params()
        self.table.update_table(self.ebp)


    def update_lattice(self):

        self.mw.lattice = deepcopy(self.mw.work_lattice)

    
    def change_value_up_down(self, values):

        for val in values:
            self.tune_elements[val[0]]['element'].__dict__[val[1]] = val[2]

        self.mw.work_lattice.lattice.update_transfer_maps()
        
        self.calc_twiss()
        self.tws_plot.update_plot(self.tws)

        self.calc_ebeam_params()
        self.table.update_table(self.ebp)
