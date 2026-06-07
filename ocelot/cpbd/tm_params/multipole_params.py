from ocelot.cpbd.tm_params.tm_params import TMParams
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.r_matrix import rot_mtx


class MultipoleParams(TMParams):
    """
    Coefficient list for the dedicated ``MultipoleTM`` kick polynomial.

    ``kn`` stores the normal multipole strengths in the order expected by the
    multipole transformation.
    """

    def __init__(self, kn) -> None:
        super().__init__()
        self.kn = kn
