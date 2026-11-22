from ocelot.cpbd.tm_params.tm_params import TMParams


class MultipoleParams(TMParams):
    """
    Coefficient list for the dedicated ``MultipoleTM`` kick polynomial.

    ``kn`` stores the normal multipole strengths in the order expected by the
    multipole transformation.
    """

    def __init__(self, kn, dx=0., dy=0., tilt=0.) -> None:
        super().__init__()
        self.kn = kn
        self.dx = dx
        self.dy = dy
        self.tilt = tilt
