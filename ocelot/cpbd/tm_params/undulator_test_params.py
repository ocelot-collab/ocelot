from ocelot.cpbd.tm_params.tm_params import TMParams

class UndulatorTestParams(TMParams):
    def __init__(self, lperiod, Kx, ax) -> None:
        self.lperiod = lperiod
        self.Kx = Kx
        self.ax = ax
