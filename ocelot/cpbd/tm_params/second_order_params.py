from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.tm_utils import transfer_map_rotation, sym_matrix


class SecondOrderParams(FirstOrderParams):
    def __init__(self, R, B, T, tilt, dx, dy) -> None:
        super().__init__(R, B, tilt)
        self.T = T
        self.dx = dx
        self.dy = dy

    def get_roteted_T(self):
        return sym_matrix(transfer_map_rotation(self.R, self.T, self.tilt)[1])
