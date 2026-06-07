from ocelot.cpbd.tm_params.tm_params import TMParams


class KickParams(TMParams):
    """
    Compact strength/offset description for ``KickTM``.

    Instead of carrying a ready-made matrix, this container stores the element
    strengths and geometry offsets that ``KickTM`` uses to compute the kick
    algorithmically.
    """

    def __init__(self, dx, dy, tilt, angle, k1, k2, k3=0.):
        super().__init__()
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.dx = dx
        self.dy = dy
        self.tilt = tilt
        self.angle = angle
