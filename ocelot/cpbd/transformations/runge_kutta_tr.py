from ocelot.cpbd.transformations.runge_kutta import RungeKuttaTM


class RungeKuttaTrTM(RungeKuttaTM):
    """
    THe same method as RungeKuttaTM but only transverse dynamics is included, longitudinal dynamics is skipped
    """
    def __init__(self, s_start=0, npoints=200):
        RungeKuttaTM.__init__(self, s_start=s_start, npoints=npoints)
        self.long_dynamics = False