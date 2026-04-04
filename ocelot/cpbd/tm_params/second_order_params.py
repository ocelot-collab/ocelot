from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.tm_utils import transfer_map_rotation, sym_matrix


class SecondOrderParams(FirstOrderParams):
    """
    Parameters for a second-order map.

    ``SecondOrderParams`` extends ``FirstOrderParams`` with the quadratic
    tensor ``T`` and stores the source offsets ``dx`` / ``dy`` used when the
    atom built the offset-aware second-order map.
    """

    def __init__(self, R, B, T, tilt, dx, dy) -> None:
        super().__init__(R, B, tilt)
        self.T = T
        self.dx = dx
        self.dy = dy

    def get_rotated_T(self):
        """Return the tilted second-order tensor in symmetrized form."""
        return sym_matrix(transfer_map_rotation(self.R, self.T, self.tilt)[1])
