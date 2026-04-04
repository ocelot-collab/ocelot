from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams


class CavityParams(FirstOrderParams):
    """
    First-order cavity map plus RF metadata.

    The linear part still lives in ``R`` / ``B`` / ``tilt``. The extra fields
    carry the RF settings needed by ``CavityTM`` / ``TWCavityTM`` to compute
    energy gain and longitudinal RF terms.
    """

    def __init__(self, R, B, tilt, v, freq, phi) -> None:
        super().__init__(R, B, tilt)
        self.v = v
        self.freq = freq
        self.phi = phi
