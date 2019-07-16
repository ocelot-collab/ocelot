"""Calculation main parameters description"""

from PyQt5.QtCore import Qt

from ocelot.cpbd.beam_params import *
from ocelot.cpbd.chromaticity import *

from app.ui_forms.widgets.ebeam_table_ui import *


class EbeamWindow(QtWidgets.QWidget):

    def __init__(self, lattice, tws):
        super().__init__()

        # init user interface
        self.ui = Ui_Ebeam_Table()
        self.ui.setupUi(self)

        self.ebp_desc = []
        self.ebp_desc.append(['Energy, GeV', 3])
        self.ebp_desc.append(['Gamma', 3])
        self.ebp_desc.append(['Periods', 0])
        self.ebp_desc.append(['Tolal length, m', 3])
        self.ebp_desc.append(['Emittance, pm*rad', 3])
        self.ebp_desc.append(['Tune X / Y', 4])
        self.ebp_desc.append(['Natural chromaticity', 3])
        self.ebp_desc.append(['Alpha, 1.0e-3', 3])
        self.ebp_desc.append(['Energy spread, 1.0e-3', 3])
        self.ebp_desc.append(['Energy loss, MeV', 3])
        self.ebp_desc.append(['Jx / Jy / Je', 3])
        self.ebp_desc.append(['tau X / Y / E ms', 3])
        self.ebp_desc.append(['I1 / I2 / I3 / I4 / I5', 3])

        self.lattice = lattice
        self.tws = tws

        self.ui.table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.ui.table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.Stretch)

        self.update_table()


    def __del__(self):
        pass


    def update_table(self):

        ebp = self.calc_ebparams()

        if ebp is None:
            self.ui.table.setRowCount(1)
            self.ui.table.setSpan(0, 0, 1, 2)

            item = QtWidgets.QTableWidgetItem("No solution")
            item.setTextAlignment(Qt.AlignCenter)
            self.ui.table.setItem(0, 0, item)

            return

        n = len(ebp) 
        self.ui.table.setRowCount(n)
        self.ui.table.clearSpans()
        
        for i in range(n):
            item = QtWidgets.QTableWidgetItem(ebp[i][0])
            self.ui.table.setItem(i, 0, item)

            item = QtWidgets.QTableWidgetItem(str(ebp[i][1]))
            item.setTextAlignment(Qt.AlignCenter)
            self.ui.table.setItem(i, 1, item)


    def calc_ebparams(self):

        if self.tws.tws is None:
            return None
            
        ebparams = EbeamParams(self.lattice.lattice, self.tws.tws[0], nsuperperiod=self.lattice.nsuperperiods)

        Q = [self.tws.tws[-1].mux/2.0/np.pi*self.lattice.nsuperperiods, self.tws.tws[-1].muy/2.0/np.pi*self.lattice.nsuperperiods]
        q_ksi = natural_chromaticity(self.lattice.lattice, self.tws.tws[0], self.lattice.nsuperperiods)

        ebp = []
        ebp.append(self.format_ebp([ebparams.E], 0))
        ebp.append(self.format_ebp([ebparams.gamma], 1))
        ebp.append(self.format_ebp([self.lattice.nsuperperiods], 2))
        ebp.append(self.format_ebp([self.lattice.nsuperperiods*self.lattice.lattice.totalLen], 3))
        ebp.append(self.format_ebp([ebparams.emittance*1.0e12], 4))
        ebp.append(self.format_ebp([Q[0], Q[1]], 5))
        ebp.append(self.format_ebp([q_ksi[0], q_ksi[1]], 6))
        ebp.append(self.format_ebp([ebparams.alpha*1.0e3], 7))
        ebp.append(self.format_ebp([ebparams.sigma_e*1.0e3], 8))
        ebp.append(self.format_ebp([ebparams.U0], 9))
        ebp.append(self.format_ebp([ebparams.Jx, ebparams.Jy, ebparams.Je], 10))
        ebp.append(self.format_ebp([ebparams.tau_x, ebparams.tau_y, ebparams.tau_e], 11))
        ebp.append(self.format_ebp([ebparams.I1, ebparams.I2, ebparams.I3, ebparams.I4, ebparams.I5], 11))

        return ebp


    def format_ebp(self, params, num):

        arr = ()
        for param in params:
            arr += (str(round(param, self.ebp_desc[num][1])), )

        string = ' / '.join(arr)

        return [self.ebp_desc[num][0], string]
    