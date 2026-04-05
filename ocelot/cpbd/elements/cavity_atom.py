import numpy as np

from ocelot.cpbd.r_matrix import cavity_coupler_edge_matrix, standing_wave_cavity_matrix
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.tm_params.cavity_params import CavityParams


class CavityAtom(Element):
    """
    Standing wave RF cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    vx_{up/down}, vy_{up/down} - zero order kick of a {up/down}stream coupler
    vxx_{up/down}, vxy_{up/down} - first order kick  a {up/down}stream coupler
    """

    def __init__(self, l=0., v=0., phi=0., freq=0., vx_up=0, vy_up=0, vxx_up=0, vxy_up=0,
                 vx_down=0, vy_down=0, vxx_down=0, vxy_down=0, eid=None, **kwargs):
        kwargs.setdefault('width', 0.1)       # Bends are usually wider than 0.05
        kwargs.setdefault('height', 0.1)
        kwargs.setdefault('color', 'yellow') # Standard color for Dipoles
        super().__init__(eid=eid, has_edge=True, **kwargs)
        self.l = l
        self.v = v  # in GV
        self.freq = freq  # Hz
        self.phi = phi  # in grad
        self.E = 0
        self.vx_up = vx_up
        self.vy_up = vy_up
        self.vxx_up = vxx_up
        self.vxy_up = vxy_up
        self.vx_down = vx_down
        self.vy_down = vy_down
        self.vxx_down = vxx_down
        self.vxy_down = vxy_down

    def __str__(self):
        s = 'Cavity('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'v=%8.6e, ' % self.v if self.v != 0. else ""
        s += 'freq=%8.6e, ' % self.freq if np.abs(self.freq) > 1e-15 else ""
        s += 'phi=%8.6e, ' % self.phi if np.abs(self.phi) > 1e-15 else ""
        s += "vx_up=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vx_up) if self.vx_up != 0. else ""
        s += "vy_up=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vy_up) if self.vy_up != 0. else ""
        s += "vxx_up=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vxx_up) if self.vxx_up != 0. else ""
        s += "vxy_up=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vxy_up) if self.vxy_up != 0. else ""
        s += "vx_down=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vx_down) if self.vx_down != 0. else ""
        s += "vy_down=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vy_down) if self.vy_down != 0. else ""
        s += "vxx_down=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vxx_down) if self.vxx_down != 0. else ""
        s += "vxy_down=({num.real:8.6e} + {num.imag:8.6e}j), ".format(num=self.vxy_down) if self.vxy_down != 0. else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s

    def _R_edge_matrix(self, energy: float, vxx: float, vxy: float):
        return cavity_coupler_edge_matrix(v=self.v, phi=self.phi, vxx=vxx, vxy=vxy, energy=energy, xp=np)

    def _R_main_matrix(self, energy: float, length: float):
        voltage = self.v * length / self.l if self.l != 0 else self.v
        return standing_wave_cavity_matrix(length, voltage=voltage, energy=energy, freq=self.freq, phi=self.phi, xp=np)

    def linear_r_main(self, energy: float = 0.0, delta_length: float = None, xp=np):
        length = self.l if delta_length is None else delta_length
        voltage = self.v * length / self.l if self.l != 0 else self.v
        return standing_wave_cavity_matrix(length, voltage=voltage, energy=energy, freq=self.freq, phi=self.phi, xp=xp)

    def linear_r_entrance(self, energy: float = 0.0, xp=np):
        return cavity_coupler_edge_matrix(v=self.v, phi=self.phi, vxx=self.vxx_up, vxy=self.vxy_up, energy=energy, xp=xp)

    def linear_r_exit(self, energy: float = 0.0, xp=np):
        return cavity_coupler_edge_matrix(v=self.v, phi=self.phi, vxx=self.vxx_down, vxy=self.vxy_down, energy=energy, xp=xp)

    def kick_b(self, v, vx, vy, phi, energy):
        phi = phi * np.pi / 180.
        dxp = (vx * v * np.exp(1j * phi)).real / energy
        dyp = (vy * v * np.exp(1j * phi)).real / energy
        b = np.array([[0.], [dxp], [0.], [dyp], [0.], [0.]])
        return b

    def create_first_order_main_params(self, energy: float, delta_length: float) -> FirstOrderParams:
        R = self.linear_r_main(energy=energy, delta_length=delta_length, xp=np)
        B = self._default_B(R)
        return FirstOrderParams(R, B, self.tilt)

    def create_first_order_entrance_params(self, energy: float, delta_length: float) -> FirstOrderParams:
        R = self.linear_r_entrance(energy=energy, xp=np)
        B = self.kick_b(self.v, self.vx_up, self.vy_up, self.phi, energy)
        return FirstOrderParams(R, B, self.tilt)

    def create_first_order_exit_params(self, energy: float, delta_length: float) -> FirstOrderParams:
        R = self.linear_r_exit(energy=energy, xp=np)
        B = self.kick_b(self.v, self.vx_down, self.vy_down, self.phi, energy)
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
