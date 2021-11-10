from ocelot.cpbd.tm_params.tm_params import TMParams


class RungeKuttaParams(TMParams):
    def __init__(self, mag_field) -> None:
        super().__init__()
        self.mag_field = mag_field
    