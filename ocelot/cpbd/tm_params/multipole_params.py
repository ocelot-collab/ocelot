from ocelot.cpbd.tm_params.tm_params import TMParams
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.r_matrix import rot_mtx


class MultipoleParams(TMParams):
    def __init__(self, kn) -> None:
        super().__init__()
        self.kn = kn
