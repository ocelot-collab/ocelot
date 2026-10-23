import logging

import numpy as np

from ocelot.cpbd.transformations.tw_cavity import TWCavityTM
from ocelot.cpbd.rf_maps import rf_long_r_terms
from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.tm_params.cavity_params import CavityParams

logger = logging.getLogger(__name__)


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

    def _R_main_matrix(self, energy: float, length: float):

        def tw_cavity_R_z(z, V, E, freq, phi=0.):
            """
            Body matrix of a travelling-wave structure.

            :param z: length of the (slice of the) structure [m]
            :param V: voltage of the (slice of the) structure [GV]
            :param E: initial energy [GeV]
            :param freq: frequency [Hz]
            :param phi: phase [deg]
            :return: matrix
            """
            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            # log1p keeps the full precision of the small-gradient limit
            r12 = z * E / de * np.log1p(de / E) if de != 0 else z
            r22 = E / (E + de)
            r55_cor, r56, r65, r66 = rf_long_r_terms(z=z, V=V, E=E, freq=freq, phi=phi)
            cav_matrix = np.array([[1, r12, 0., 0., 0., 0.],
                                   [0, r22, 0., 0., 0., 0.],
                                   [0., 0., 1, r12, 0., 0.],
                                   [0., 0., 0, r22, 0., 0.],
                                   [0., 0., 0., 0., 1. + r55_cor, r56],
                                   [0., 0., 0., 0., r65, r66]]).real
            return cav_matrix

        if self.v == 0.:
            R = uni_matrix(length, 0., hx=0., sum_tilts=self.tilt, energy=energy)
        else:
            R = tw_cavity_R_z(length, V=self.v * length / self.l, E=energy, freq=self.freq,
                           phi=self.phi)
        return R

    def _R_entrance(self, z, V, E, phi=0.):
        """
        Entrance edge lens, -gamma'/(2 gamma0).

        `E` is the reference energy in front of the structure.
        """
        phi = phi * np.pi / 180.
        de = V * np.cos(phi)
        r = np.eye(6)
        r[1, 0] = -de / z / 2. / E
        r[3, 2] = r[1, 0]
        return r

    def _R_exit(self, z, V, E, phi=0.):
        """
        Exit edge lens, +gamma'/(2 gamma1).

        `E` is the reference energy *behind* the structure: the transfer-map
        chain has already added the energy gain of the body before this edge
        map is built, so it must not be added a second time here.
        """
        phi = phi * np.pi / 180.
        de = V * np.cos(phi)
        r = np.eye(6)
        r[1, 0] = +de / z / 2. / E
        r[3, 2] = r[1, 0]
        return r

    def create_first_order_main_params(self, energy: float, delta_length: float) -> FirstOrderParams:
        R = self._R_main_matrix(energy=energy, length=delta_length if delta_length is not None else self.l)
        B = self._default_B(R)
        return FirstOrderParams(R, B, self.tilt)

    def create_first_order_entrance_params(self, energy: float, delta_length: float) -> FirstOrderParams:
        length = delta_length if delta_length is not None else self.l
        R = self._R_entrance(z=length, V=self.v * length / self.l, E=energy, phi=self.phi)
        B = self._default_B(R)
        return FirstOrderParams(R, B, self.tilt)

    def create_first_order_exit_params(self, energy: float, delta_length: float) -> FirstOrderParams:
        length = delta_length if delta_length is not None else self.l
        R = self._R_exit(z=length, V=self.v * length / self.l, E=energy, phi=self.phi)
        B = self._default_B(R)
        return FirstOrderParams(R, B, self.tilt)

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
