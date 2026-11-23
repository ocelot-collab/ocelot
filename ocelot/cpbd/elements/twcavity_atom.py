import numpy as np

from ocelot.common.globals import m_e_GeV
from ocelot.cpbd.rf_maps import rf_long_r_terms
from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.tm_params.cavity_params import CavityParams

class TWCavityAtom(Element):
    """
    Traveling wave cavity (phase velocity equal to the speed of light).

    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]

    The body of a travelling-wave structure has no ponderomotive focusing: it
    only damps adiabatically, and the net radial focusing comes from the
    entrance/exit edge lenses -gamma'/(2 gamma0) and +gamma'/(2 gamma1)
    (Rosenzweig & Serafini, PRE 49, 1599 (1994)). The longitudinal map is the
    same as for a standing wave and is shared through ocelot.cpbd.rf_maps.
    """
    def __init__(self, l=0., v=0., phi=0., freq=0., eid=None, **kwargs):
        super().__init__(eid=eid, has_edge=True, **kwargs)
        self.l = l
        self.v = v  # in GV
        self.freq = freq  # Hz
        self.phi = phi  # in grad
        self.E = 0

    def __str__(self):
        s = 'TWCavity('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'v=%8.6e, ' % self.v if self.v != 0. else ""
        s += 'freq=%8.6e, ' % self.freq if np.abs(self.freq) > 1e-15 else ""
        s += 'phi=%8.6e, ' % self.phi if np.abs(self.phi) > 1e-15 else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s

    def _R_main_matrix(self, energy: float, length: float) -> np.ndarray:
        """Transport with constant transverse momentum through a linear energy ramp."""
        if self.v == 0.:
            return uni_matrix(length, 0., hx=0., sum_tilts=self.tilt, energy=energy)
        if self.l == 0:
            raise ValueError("A powered travelling-wave cavity requires nonzero length")
        phi = np.deg2rad(self.phi)
        voltage = self.v * length / self.l
        r55_cor, r56, r65, r66 = rf_long_r_terms(length, voltage, energy, self.freq, phi)

        g0 = energy / m_e_GeV
        dgamma = voltage * np.cos(phi) / m_e_GeV
        g1 = g0 + dgamma
        p0 = np.sqrt((g0 - 1.) * (g0 + 1.))
        p1 = np.sqrt((g1 - 1.) * (g1 + 1.))
        # R12 = p0 integral(ds/p). The integral is a difference of rapidities
        # divided by dgamma. Rationalize p1-p0 for very small energy gains.
        if dgamma:
            dp = dgamma * (g1 + g0) / (p1 + p0)
            r12 = length * p0 * np.log1p((dgamma + dp) / (g0 + p0)) / dgamma
        else:
            r12 = length
        matrix = np.eye(6)
        matrix[0, 1] = matrix[2, 3] = r12
        matrix[1, 1] = matrix[3, 3] = p0 / p1
        matrix[4:, 4:] = [[1. + r55_cor, r56], [r65, r66]]
        return matrix

    def _edge_params(self, energy: float, sign: float) -> FirstOrderParams:
        """Thin edge lens at the physical face, using the energy at that face."""
        matrix = np.eye(6)
        if self.v:
            if self.l == 0 or energy <= m_e_GeV:
                raise ValueError("A powered travelling-wave cavity requires nonzero length and energy above rest energy")
            gradient = self.v * np.cos(np.deg2rad(self.phi)) / self.l
            matrix[1, 0] = matrix[3, 2] = sign * gradient / (2. * energy)
        return FirstOrderParams(matrix, self._default_B(matrix), self.tilt)

    def create_first_order_main_params(self, energy: float, delta_length: float | None) -> FirstOrderParams:
        matrix = self._R_main_matrix(energy, self.l if delta_length is None else delta_length)
        return FirstOrderParams(matrix, self._default_B(matrix), self.tilt)

    def create_first_order_entrance_params(self, energy: float, delta_length: float | None) -> FirstOrderParams:
        return self._edge_params(energy, -1.)

    def create_first_order_exit_params(self, energy: float, delta_length: float | None) -> FirstOrderParams:
        return self._edge_params(energy, +1.)

    def create_cavity_tm_main_params(self, energy: float, delta_length: float) -> CavityParams:
        fo_params = self.create_first_order_main_params(energy, delta_length)
        return CavityParams(R=fo_params.R, B=fo_params.B, tilt=self.tilt, v=self.v, freq=self.freq, phi=self.phi)

    def create_cavity_tm_entrance_params(self, energy: float, delta_length: float) -> CavityParams:
        fo_params = self.create_first_order_entrance_params(energy, delta_length)
        return CavityParams(R=fo_params.R, B=fo_params.B, tilt=self.tilt, v=self.v, freq=self.freq, phi=self.phi)

    def create_cavity_tm_exit_params(self, energy: float, delta_length: float) -> CavityParams:
        fo_params = self.create_first_order_exit_params(energy, delta_length)
        return CavityParams(R=fo_params.R, B=fo_params.B, tilt=self.tilt, v=self.v, freq=self.freq, phi=self.phi)

    def create_delta_e(self, total_length, delta_length=None) -> float:
        if delta_length is not None:
            phase_term = np.cos(self.phi * np.pi / 180.)
            return phase_term * self.v * delta_length / total_length if total_length != 0 else phase_term * self.v
        else:
            return self.v * np.cos(self.phi * np.pi / 180.)
