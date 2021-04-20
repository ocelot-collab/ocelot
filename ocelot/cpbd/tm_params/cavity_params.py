from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams

class CavityParams(FirstOrderParams):
    def __init__(self, R, B, tilt, v, freq, phi) -> None:
        super().__init__(R, B, tilt)
        self.v = v
        self.freq = freq
        self.phi = phi
