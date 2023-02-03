from ocelot.cpbd.transformations.runge_kutta import RungeKuttaTM, TMTypes


class RungeKuttaTrTM(RungeKuttaTM):
    """
    THe same method as RungeKuttaTM but only transverse dynamics is included, longitudinal dynamics is skipped
    """
    def __init__(self, create_tm_param_func, delta_e_func, tm_type: TMTypes, length: float, delta_length: float = 0.0, **params) -> None:
        super().__init__(create_tm_param_func, delta_e_func, tm_type, length, delta_length, **params)
        self.long_dynamics = False