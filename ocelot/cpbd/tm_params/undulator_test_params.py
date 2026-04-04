from ocelot.cpbd.tm_params.tm_params import TMParams

class UndulatorTestParams(TMParams):
    """
    Minimal undulator metadata used by the legacy/test undulator map.

    This is not a general field container; it only carries the quantities that
    ``UndulatorTestTM`` integrates with its simplified test algorithm.
    """

    def __init__(self, lperiod, Kx, ax) -> None:
        self.lperiod = lperiod
        self.Kx = Kx
        self.ax = ax
