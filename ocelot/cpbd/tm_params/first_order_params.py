import numpy as np
from ocelot.cpbd.tm_params.tm_params import TMParams
from ocelot.cpbd.r_matrix import rot_mtx


class FirstOrderParams(TMParams):
    def __init__(self, R, B, tilt) -> None:
        super().__init__()
        self.R = R
        self.B = B
        self.tilt = tilt

    def get_rotated_R(self):
        return np.dot(np.dot(rot_mtx(-self.tilt), self.R), rot_mtx(self.tilt))

    def multiply(self, tm_params: 'FirstOrderParams'):
        R = np.dot(self.R, tm_params.R)
        B = np.dot(self.R, tm_params.B) + self.B
        return FirstOrderParams(R, B, self.tilt)