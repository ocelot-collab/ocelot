import numpy as np

from ocelot.cpbd.transformations.transfer_map import TransferMap, TMTypes
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.rf_maps import rf_longitudinal_kick


class CavityTM(TransferMap):
    """Track through RF edge lenses and a cavity body.

    The element atom supplies the linear transverse and time-of-flight maps.
    Inside the body both cavity families use the particle's entrance phase
    for its energy gain and a second-order correction to arrival time.
    """

    def __init__(self, create_tm_param_func, delta_e_func, tm_type: TMTypes, length: float,
                 delta_length: float) -> None:
        super().__init__(create_tm_param_func, delta_e_func, tm_type, length, delta_length=delta_length)

    @classmethod
    def from_element(cls, element: Element, tm_type: TMTypes = TMTypes.MAIN, delta_l=None, **params):
        return cls.create(entrance_tm_params_func=element.create_cavity_tm_entrance_params,
                          delta_e_func=element.create_delta_e,
                          main_tm_params_func=element.create_cavity_tm_main_params,
                          exit_tm_params_func=element.create_cavity_tm_exit_params,
                          tm_type=tm_type, length=element.l, delta_length=delta_l)

    def map4cav(self, X: np.ndarray, E: float, delta_length: float | None, length: float) -> np.ndarray:
        params = self.get_params(E)
        z = length if delta_length is None else delta_length
        voltage = params.v * z / length if length else params.v
        tau, delta = X[4].copy(), X[5].copy()
        X = self.mul_p_array(X, energy=E)
        return rf_longitudinal_kick(X, tau, delta, E, voltage, z, params.freq, np.deg2rad(params.phi))

    def map_function(self, X, energy: float):
        if self.tm_type == TMTypes.MAIN:
            return self.map4cav(X, energy, self.delta_length, self.length)
        else:
            return self.mul_p_array(X, energy=energy)
