import numpy as np

from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.elements.element import Element


class Matrix(Element):
    """
    Matrix element

    l = 0 - m, length of the matrix element
    r = np.zeros((6, 6)) - R - elements, first order
    t = np.zeros((6, 6, 6)) - T - elements, second order
    delta_e = 0 - GeV, energy gain along the matrix element
    """

    def __init__(self, l=0., delta_e=0, eid=None, **kwargs):
        Element.__init__(self, eid)
        self.l = l

        self.r = np.zeros((6, 6))
        self.t = np.zeros((6, 6, 6))
        # zero order elements - test mode, not implemented yet
        self.b = np.zeros((6, 1))

        for y in kwargs:
            # decode first order arguments in format RXX or rXX where X is number from 1 to 6
            if "r" in y[0].lower() and len(y) > 2:
                if "m" in y[1].lower() and len(y) == 4 and y[2:].isdigit() and (11 <= int(y[2:]) <= 66):
                    self.r[int(y[2]) - 1, int(y[3]) - 1] = float(kwargs[y])
                if len(y) == 3 and y[1:].isdigit() and (11 <= int(y[1:]) <= 66):
                    self.r[int(y[1]) - 1, int(y[2]) - 1] = float(kwargs[y])

            # decode second order arguments in format TXXX or tXXX where X is number from 1 to 6
            if "t" in y[0].lower() and len(y) == 4 and y[1:].isdigit() and (111 <= int(y[1:]) <= 666):
                self.t[int(y[1]) - 1, int(y[2]) - 1, int(y[3]) - 1] = float(kwargs[y])

            # decode zero order arguments in format BX or bX where X is number from 1 to 6
            if "b" in y[0].lower() and len(y) == 2 and y[1:].isdigit() and (1 <= int(y[1:]) <= 6):
                self.b[int(y[1]) - 1, 0] = float(kwargs[y])
        self.delta_e = delta_e

    def __str__(self):
        s = 'Matrix : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l =%8.5f m\n' % self.l
        s += 'R = \n'
        for i in range(6):
            for j in range(6):
                s += '%11.6f' % (self.r[i, j])
            s += "\n"
        return s

    def create_r_matrix(self):
        rm = np.eye(6)
        rm = self.r

        def r_matrix(z, l, rm):
            if z < l:
                r_z = uni_matrix(z, 0, hx=0)
            else:
                r_z = rm
            return r_z

        r_z_e = lambda z, energy: r_matrix(z, self.l, rm)
        return r_z_e

    def get_T_z_e_func(self):
        return lambda z, energy: self.t

    def _set_general_tm_parameter(self):
        self.transfer_map.delta_e = self.delta_e
        self.transfer_map.B_z = lambda z, energy: self.b
        self.transfer_map.B = lambda energy: self.b
        super()._set_general_tm_parameter()