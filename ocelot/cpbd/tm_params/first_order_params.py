import numpy as np
from ocelot.cpbd.tm_params.tm_params import TMParams
from ocelot.cpbd.r_matrix import rot_mtx


class FirstOrderParams(TMParams):
    """
    Parameters for a linear first-order map.

    Attributes
    ----------
    R
        6x6 transfer matrix before tilt rotation is applied.
    B
        6x1 additive offset term for the linear map.
    tilt
        Roll angle used when the transformation requests a rotated matrix.
    """

    def __init__(self, R, B, tilt) -> None:
        super().__init__()
        self.R = R
        self.B = B
        self.tilt = tilt
        self._rotated_R = None

    def get_rotated_R(self):
        """Return ``R`` in the tilted element frame expected by the map."""
        if self._rotated_R is None:
            if self.tilt == 0:
                self._rotated_R = self.R
            else:
                self._rotated_R = np.dot(np.dot(rot_mtx(-self.tilt), self.R), rot_mtx(self.tilt))
        return self._rotated_R.copy()

    def multiply(self, tm_params: 'FirstOrderParams'):
        """Compose two linear parameter sets into one equivalent map."""
        R = np.dot(self.R, tm_params.R)
        B = np.dot(self.R, tm_params.B) + self.B
        return FirstOrderParams(R, B, self.tilt)
