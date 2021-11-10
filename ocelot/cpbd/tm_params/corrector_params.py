from ocelot.cpbd.tm_params.second_order_params import SecondOrderParams


class CorrectorSecondOrderParams(SecondOrderParams):
    def __init__(self, R, B, T, tilt, dx, dy, angle) -> None:
        super().__init__(R, B, T, tilt, dx, dy)
        self.angle = angle
