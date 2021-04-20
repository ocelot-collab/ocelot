import numpy as np

from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.second_order_params import SecondOrderParams
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams


class MatrixAtom(Element):
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

    def create_delta_e(self, total_length, delta_length=0.0):
        return self.delta_e

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

    def create_first_order_main_params(self, energy: float, delta_length: float) -> FirstOrderParams:
        rm = self.r
        if delta_length != None and delta_length < self.l:
            R = uni_matrix(delta_length, 0, hx=0)
        else:
            R = rm
        return FirstOrderParams(R=R, B=self.b, tilt=self.tilt)

    def create_second_order_main_params(self, energy: float, delta_length: float) -> SecondOrderParams:
        fo_params = self.create_first_order_main_params(energy, delta_length=delta_length)
        return SecondOrderParams(R=fo_params.R, B=fo_params.B, tilt=fo_params.tilt, dx=self.dx, dy=self.dy, T=self.t)
