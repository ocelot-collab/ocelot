import numpy as np

from ocelot.cpbd.transformations.transfer_map import TransferMap, TMTypes
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.rf_maps import rf_longitudinal_kick


class TWCavityTM(TransferMap):
    """[summary]
    Implementation of the Traveling Wave Cavity Transformation.
    The concrete element atom have to implement:
    create_cavity_tm_main_params(self, energy: float, delta_length: float) -> CavityParams
    create_cavity_tm_entrance_params(self, energy: float, delta_length: float) -> CavityParams
    create_cavity_tm_exit_params(self, energy: float, delta_length: float) -> CavityParams

    The structure differs from a standing-wave cavity only in its transverse
    block, which the atom supplies through the R matrix. The longitudinal map
    is identical for both families and lives in ocelot.cpbd.rf_maps.
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

    def map4cav(self, X, E, delta_length, length):
        params = self.get_params(E)

        # Effective cavity voltage and length slice
        if delta_length is not None:
            V = params.v * delta_length / length if length != 0 else params.v
            z = delta_length
        else:
            V = params.v
            z = length

        phi = params.phi * np.pi / 180.0

        # Save initial longitudinal coordinates
        X4 = np.copy(X[4])  # sigma
        X5 = np.copy(X[5])  # delta

        # Apply the linear part (transverse damping + longitudinal block)
        X = self.mul_p_array(X, energy=E)

        # Longitudinal RF map, shared with the standing-wave cavity
        X = rf_longitudinal_kick(X, X4, X5, E=E, V=V, z=z, freq=params.freq, phi=phi)

        return X

    def map_function(self, X, energy: float):
        if self.tm_type == TMTypes.MAIN:
            return self.map4cav(X, energy, self.delta_length, self.length)
        else:
            return self.mul_p_array(X, energy=energy)
