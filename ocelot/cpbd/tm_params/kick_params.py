from ocelot.cpbd.tm_params.tm_params import TMParams


class KickParams(TMParams):
    def __init__(self, dx, dy, tilt, angle, k1, k2, k3=0.):
        super().__init__()
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.dx = dx
        self.dy = dy
        self.tilt = tilt
        self.angle = angle