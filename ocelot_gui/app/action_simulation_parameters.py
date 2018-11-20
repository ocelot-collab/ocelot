from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QWidget

from copy import deepcopy

from app.ui_forms.widgets.ebeam_table import *

from ocelot.cpbd.e_beam_params import *
from ocelot.cpbd.chromaticity import *


class CalcParameters(QWidget):

    def __init__(self, MainWindow):
        super().__init__()
        
        self.setWindowTitle('Main parameters')
        self.resize(500, 500)
        
        self.mw = MainWindow

        self.init_form()

    
    def init_form(self):
        
        self.mw.lattice.init_lattice()

        self.table = EbeamTable()

        self.calc_twiss()
        self.calc_ebeam_params()
        self.table.update_table(self.ebp)

        layout = QtGui.QGridLayout()
        layout.addWidget(self.table.table)

        self.setLayout(layout)


    def calc_twiss(self):

        if self.mw.lattice.periodic_solution:
            tws0 = periodic_twiss(Twiss(), lattice_transfer_map(self.mw.lattice.lattice, self.mw.lattice.tws0.E))
        else:
            tws0 = self.mw.lattice.tws0

        if tws0 is not None:
            self.tws = twiss(self.mw.lattice.lattice, tws0=tws0, nPoints=None)
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
        self.ebp.append(['I1 / I2 / I3 / I4 / I5', '', 3])


    def calc_ebeam_params(self):

        self.init_ebp()

        if self.tws is None:
            self.ebp[0] = None
            return
        
        (I1, I2, I3, I4, I5) = radiation_integrals(self.mw.lattice.lattice, self.tws[0], self.mw.lattice.nsuperperiods)
        if I2 == 0.0:
            self.ebp[0] = None
            return

        Q = [self.tws[-1].mux/2.0/np.pi*self.mw.lattice.nsuperperiods, self.tws[-1].muy/2.0/np.pi*self.mw.lattice.nsuperperiods]
        q_ksi = natural_chromaticity(self.mw.lattice.lattice, self.tws[0], self.mw.lattice.nsuperperiods)
        alpha = I1 / (self.mw.lattice.nsuperperiods * self.mw.lattice.lattice.totalLen) * 1.0e3
        Jx, Jy, Je = 1.0 - I4 / I2, 1.0, 2.0 + I4 / I2

        self.ebp[3][1] = self.mw.lattice.nsuperperiods
        self.ebp[4][1] = round(self.mw.lattice.nsuperperiods*self.mw.lattice.lattice.totalLen, self.ebp[4][2])
        self.ebp[6][1] = str(round(Q[0], self.ebp[6][2])) + ' / ' + str(round(Q[1], self.ebp[6][2]))
        self.ebp[7][1] = str(round(q_ksi[0], self.ebp[7][2])) + ' / ' + str(round(q_ksi[1], self.ebp[7][2]))
        self.ebp[8][1] = round(alpha, self.ebp[8][2])
        self.ebp[11][1] = str(round(Jx, self.ebp[11][2])) + ' / ' + str(round(Jy, self.ebp[11][2])) + ' / ' + str(round(Je, self.ebp[11][2]))

        if self.mw.lattice.tws0.E > 0.0:
            
            gamma = self.mw.lattice.tws0.E/m_e_GeV
            emittance = Cq * gamma**2 * I5 / Jx / I2 * 1.0e12
            u0 = Cgamma * (self.mw.lattice.tws0.E * 1000.0)**4 * I2 / 2.0 / np.pi
            sigma_e = gamma * np.sqrt(Cq * I3 / Je / I2) * 1.0e3
            tau0 = 2.0 * self.mw.lattice.tws0.E * 1000.0 * self.mw.lattice.nsuperperiods * self.mw.lattice.lattice.totalLen / speed_of_light / u0 * 1.0e3

            self.ebp[1][1] = round(self.mw.lattice.tws0.E, self.ebp[1][2])
            self.ebp[2][1] = round(gamma, self.ebp[12][2])
            self.ebp[5][1] = round(emittance, self.ebp[5][2])
            self.ebp[9][1] = round(sigma_e, self.ebp[9][2])
            self.ebp[10][1] = round(u0, self.ebp[10][2])
            self.ebp[12][1] = str(round(tau0/Jx, self.ebp[12][2])) + ' / ' + str(round(tau0/Jy, self.ebp[12][2])) + ' / ' + str(round(tau0/Je, self.ebp[12][2]))
            
        self.ebp[13][1] = str(round(I1, self.ebp[13][2])) + ' / ' + str(round(I2, self.ebp[13][2])) + ' / ' + str(round(I3, self.ebp[13][2])) + ' / ' + str(round(I4, self.ebp[13][2])) + ' / ' + str(round(I5, self.ebp[13][2]))
