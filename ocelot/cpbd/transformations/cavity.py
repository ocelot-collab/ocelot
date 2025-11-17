import numpy as np

from ocelot.cpbd.transformations.transfer_map import TransferMap, TMTypes
from ocelot.cpbd.elements.element import Element
from ocelot.common.globals import speed_of_light, m_e_GeV


class CavityTM(TransferMap):
    """[summary]
    Implementation of the Cavity Transformation.
    The concrete element atom have to implement: 
    create_cavity_tm_main_params(self, energy: float, delta_length: float) -> CavityParams
    create_cavity_tm_entrance_params(self, energy: float, delta_length: float) -> CavityParams
    create_cavity_tm_exit_params(self, energy: float, delta_length: float) -> CavityParams
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

        # Initial gamma/beta
        beta0 = 1.0
        igamma2 = 0.0
        g0 = 1e10  # large gamma if E == 0
        if E != 0.0:
            g0 = E / m_e_GeV
            igamma2 = 1.0 / (g0 * g0)
            beta0 = np.sqrt(1.0 - igamma2)

        phi = params.phi * np.pi / 180.0

        # Save initial longitudinal coordinates
        X4 = np.copy(X[4])  # sigma
        X5 = np.copy(X[5])  # delta

        # Apply transverse/linear RF focusing (R, T already implemented elsewhere)
        X = self.mul_p_array(X, energy=E)

        # Energy gain of reference particle
        delta_e = V * np.cos(phi)
        E1 = E + delta_e

        # Default: drift-like second order (T566^drift, T555=T556=0)
        T566 = 1.5 * z * igamma2 / (beta0 ** 3)  # = 3L/(2 beta^3 gamma^2)
        T556 = 0.0
        T555 = 0.0

        # If final energy would be non-physical, just return drift map
        if E1 <= 0.0:
            X[4] += T566 * X5 * X5  # T556=T555=0 here
            return X

        # RF wavenumber
        k = 2.0 * np.pi * params.freq / speed_of_light

        # Final gamma/beta
        g1 = E1 / m_e_GeV
        beta1 = np.sqrt(1.0 - 1.0 / (g1 * g1))

        # Longitudinal energy coordinate (delta) update:
        # δ1 = E0 β0 / (E1 β1) δ0 + V β0 / (E1 β1) [cos(φ0 - k β0 σ0) - cos φ0]
        X[5] = (
                X5 * E * beta0 / (E1 * beta1)
                + V * beta0 / (E1 * beta1) * (np.cos(-X4 * beta0 * k + phi) - np.cos(phi))
        )

        # General RF second-order coefficients
        dgamma = V / m_e_GeV
        dg = g1 - g0

        # Threshold for "small energy gain" (relative)
        eps_rel = 1e-8

        if abs(dg) < eps_rel * abs(g0):
            # --- Small Δγ limit (γ1 ≈ γ0) ---
            # T566 tends to drift value (already set above).
            if abs(np.cos(phi)) < 1e-3:
                # Zero-crossing cavity: use analytic limits for T555 and T556
                # T556_lim = 3/2 * L * k * dgamma / (γ0^3 β0^3)
                T556 = 1.5 * z * k * dgamma / (beta0 ** 3 * g0 ** 3)

                # T555_lim = 1/2 * L * k^2 * dgamma^2 / (γ0^4 β0^3)
                T555 = 0.5 * z * k * k * dgamma * dgamma / (beta0 ** 3 * g0 ** 4)
            else:
                # Small-gradient limit away from zero-crossing:
                # cavity behaves like a drift in second order
                T556 = 0.0
                T555 = 0.0
        else:
            # --- General case: full RF formulas, valid for any sign of ΔE ---
            T566 = ( z * (beta0 ** 3 * g0 ** 3 - beta1 ** 3 * g1 ** 3)
                    / (2.0 * beta0 * beta1 ** 3 * g0 * (g0 - g1) * g1 ** 3))

            T556 = (beta0 * k * z * dgamma * g0  * (beta1 ** 3 * g1 ** 3 + beta0 * (g0 - g1 ** 3)) * np.sin(phi) /
                    (beta1 ** 3 * g1 ** 3 * (g0 - g1) ** 2))

            T555 = ( beta0 ** 2 * k ** 2 * z * dgamma / 2.0
                    * ( dgamma * ( 2.0 * g0 * g1 ** 3 * (beta0 * beta1 ** 3 - 1.0)
                                    + g0 ** 2 + 3.0 * g1 ** 2 - 2.0 ) / (beta1 ** 3 * g1 ** 3 * (g0 - g1) ** 3)  * np.sin(phi) ** 2
                            - (g1 * g0 * (beta1 * beta0 - 1.0) + 1.0)
                            / (beta1 * g1 * (g0 - g1) ** 2)
                            * np.cos(phi)
                    )
            )

        # Apply second-order longitudinal map (σ update)
        X[4] += T566 * X5 * X5 + T556 * X4 * X5 + T555 * X4 * X4

        return X

    def map_function(self, X, energy: float):
        if self.tm_type == TMTypes.MAIN:
            return self.map4cav(X, energy, self.delta_length, self.length)
        else:
            return self.mul_p_array(X, energy=energy)
