# Explicit re-exports used by existing simulation scripts.
from ocelot.cpbd.elements import Quadrupole as Quadrupole

from ocelot.cpbd.beam import Twiss
from ocelot.cpbd.match import closed_orbit
from ocelot.cpbd.elements import Cavity, Hcor, TWCavity
from ocelot.cpbd.beam import Particle
from ocelot.common.globals import m_e_GeV
from ocelot.cpbd.transformations.linearization import map_derivatives
import copy
import os
import numpy as np
import json
import pandas as pd

import logging

logger = logging.getLogger(__name__)


class MeasureResponseMatrix:
    def __init__(self, lattice, hcors, vcors, bpms):
        self.lat = lattice
        self.hcors = hcors
        self.vcors = vcors
        self.bpms = bpms

    def calculate(self):
        pass

    def read_virtual_orbit(self, p_init: Particle | None = None,
                           write2bpms: bool = True) -> tuple[np.ndarray, np.ndarray]:
        """Return horizontal and vertical positions [m] in the requested BPM order.

        A lattice BPM reads at its centre. Lattice order distinguishes BPMs
        separated by a thin element at the same s. A detached position probe
        uses its s [m] from the lattice entrance, downstream of thin elements
        at that position. With no launch particle, use the closed orbit.
        """
        self.particle0 = closed_orbit(self.lat) if p_init is None else p_init
        p = copy.copy(self.particle0)
        sequence = self.lat.sequence
        boundaries = np.r_[0., np.cumsum([elem.l for elem in sequence])]
        samples = []
        for row, bpm in enumerate(self.bpms):
            if self.lat.find_element_indices(bpm) or not hasattr(bpm, "s"):
                index = self.lat.resolve_element_index(bpm)
                offset = bpm.l / 2.
            else:
                if not 0. <= bpm.s <= boundaries[-1]:
                    raise ValueError("BPM position must lie inside the lattice")
                index = min(np.searchsorted(boundaries, bpm.s, side="right") - 1, len(sequence))
                offset = bpm.s - boundaries[index]
            samples.append((index, offset, row))

        X, Y = np.empty(len(samples)), np.empty(len(samples))
        current = 0
        for index, target, row in sorted(samples):
            while current < index:
                for tm in sequence[current].tms:
                    tm.apply(p)
                current += 1
            reading = p
            if target > 0.:
                # A measurement inside a magnet must not change its downstream
                # trajectory. Splitting a numerical map can change its kicks.
                reading = copy.copy(p)
                for tm in sequence[index].get_section_tms(target):
                    tm.apply(reading)
            if write2bpms:
                bpm = self.bpms[row]
                bpm.x, bpm.y, bpm.E, bpm.p = reading.x, reading.y, reading.E, reading.p
            X[row], Y[row] = reading.x, reading.y
        return X, Y

    def read_virtual_dispersion(self, E0: float, write2bpms: bool = True) -> tuple[np.ndarray, np.ndarray]:
        """Track BPM dispersion using entrance Particle.p values of 0 and 0.01.

        Both particles use the same entrance reference energy E0 [GeV].
        The denominator is the entrance offset, including in accelerating lines.
        """
        dp = 0.01
        X0, Y0 = self.read_virtual_orbit(p_init=Particle(E=E0), write2bpms=write2bpms)
        X1, Y1 = self.read_virtual_orbit(p_init=Particle(p=dp, E=E0), write2bpms=False)
        Dx = (X1 - X0) / dp
        Dy = (Y1 - Y0) / dp
        if write2bpms:
            for bpm, dx, dy in zip(self.bpms, Dx, Dy):
                bpm.Dx = dx
                bpm.Dy = dy
        return Dx, Dy

    def _calculate_kick_response(self, read_bpms, angle_step):
        """Forward differences around the current settings. Only baseline readings are saved."""
        correctors = [*self.hcors, *self.vcors]
        self._device_indices(correctors)
        baseline = np.concatenate(read_bpms(write2bpms=True))
        response = np.empty((baseline.size, len(correctors)))
        for j, cor in enumerate(correctors):
            original_angle = cor.angle
            try:
                cor.angle = original_angle + angle_step
                measured = np.concatenate(read_bpms(write2bpms=False))
                response[:, j] = (measured - baseline) / angle_step
            finally:
                cor.angle = original_angle
        self.resp = response
        return response

    def _device_indices(self, devices):
        """Map unique lattice occurrence indices to their requested row or column."""
        indices = {self.lat.resolve_element_index(device): i for i, device in enumerate(devices)}
        if len(indices) != len(devices):
            raise ValueError("A response matrix device was selected more than once")
        return indices

    def _response_sections(self, first_order_only):
        """Yield occurrence index, maps, unit kick at corrector exit and BPM row at centre."""
        correctors = [*self.hcors, *self.vcors]
        columns = self._device_indices(correctors)
        rows = self._device_indices(self.bpms)
        for i, elem in enumerate(self.lat.sequence):
            row = rows.get(i)
            if row is not None and elem.l:
                half = elem.l / 2.
                yield i, elem.get_section_tms(half, first_order_only=first_order_only), None, row
                yield i, elem.get_section_tms(half, start_l=half, first_order_only=first_order_only), None, None
            else:
                kick = None
                if i in columns:
                    unit_kick = elem.element.kick_b(elem.l, elem.l, 1., elem.tilt)
                    kick = (columns[i], unit_kick[:, 0])
                maps = elem.first_order_tms if first_order_only else elem.tms
                yield i, maps, kick, row

    def _linear_response(self, energy, periodic=False):
        """Transport all unit kicks. For a ring also transport four launch coordinates."""
        ncor = len(self.hcors) + len(self.vcors)
        nbpm = len(self.bpms)
        nlaunch = 4 if periodic else 0
        sensitivities = np.zeros((6, nlaunch + ncor))
        sensitivities[:nlaunch, :nlaunch] = np.eye(nlaunch)
        readings = np.zeros((2 * nbpm, nlaunch + ncor))
        for _, maps, kick, row in self._response_sections(first_order_only=True):
            for tm in maps:
                R = tm.get_params(energy).get_rotated_R()
                sensitivities = R @ sensitivities
                energy += tm.get_delta_e()
            if kick is not None:
                column, unit_kick = kick
                sensitivities[:, nlaunch + column] += unit_kick
            if row is not None:
                readings[[row, row + nbpm]] = sensitivities[[0, 2]]

        if periodic:
            # Fixed entrance tau and delta, as in closed_orbit(). The same
            # transverse launch coordinates must recur after one turn.
            one_turn = sensitivities[:4, :4]
            if np.max(abs(np.linalg.eigvals(one_turn))) > 1. + 1.e-7:
                raise ValueError("Unstable one-turn transverse map")
            try:
                closed = np.linalg.solve(np.eye(4) - one_turn, sensitivities[:4, 4:])
            except np.linalg.LinAlgError as exc:
                raise ValueError("No unique periodic orbit response") from exc
            readings = readings[:, 4:] + readings[:, :4] @ closed
        return readings

    def _optics_coupling(self, energy):
        """Check when independent, constant-energy Courant-Snyder planes suffice."""
        longitudinal = False
        for elem in self.lat.sequence:
            for tm in elem.first_order_tms:
                R = tm.get_params(energy).get_rotated_R()
                if np.any(abs(R[:2, 2:4]) > 1.e-14) or np.any(abs(R[2:4, :2]) > 1.e-14):
                    return True, longitudinal
                gain = tm.get_delta_e()
                longitudinal = longitudinal or gain != 0. or np.any(R[:4, 4:])
                energy += gain
        return False, longitudinal

    def optical_func_params(self, tw_init=None):
        """Attach uncoupled Twiss parameters at the device centres and store the tunes."""
        optics = self.lat.periodic_twiss() if tw_init is None else tw_init
        if optics.beta_x <= 0 or optics.beta_y <= 0:
            raise ValueError("Optical functions require positive beta_x and beta_y")
        start_s, start_mux, start_muy = optics.s, optics.mux, optics.muy
        devices = self._device_indices([*self.hcors, *self.vcors, *self.bpms])
        for i, elem in enumerate(self.lat.sequence):
            if i in devices:
                centre = optics
                if elem.l:
                    for tm in elem.get_section_tms(elem.l / 2., first_order_only=True):
                        centre = tm * centre
                elem.phi_x, elem.phi_y = centre.mux, centre.muy
                elem.beta_x, elem.beta_y = centre.beta_x, centre.beta_y
                elem.E = centre.E
                elem.s = centre.s - start_s
                elem.lat_inx = i
            for tm in elem.first_order_tms:
                optics = tm * optics
        self.nu_x = (optics.mux - start_mux) / (2. * np.pi)
        self.nu_y = (optics.muy - start_muy) / (2. * np.pi)


class RingRM(MeasureResponseMatrix):
    """Linear closed-orbit response at fixed entrance tau and delta, including coupling."""

    def calculate(self, tw_init=None):
        energy = tw_init.E if tw_init is not None else 0.
        coupled, _ = self._optics_coupling(energy)
        if not coupled:
            self.optical_func_params(self.lat.periodic_twiss(tw_init))
        self.resp = self._linear_response(energy, periodic=True)
        return self.resp


class LinacOpticalRM(MeasureResponseMatrix):
    """Linear optics response, with beta and phase at devices for uncoupled lattices.

    The full optics matrices also cover coupling and RF damping. In an
    uncoupled line at constant energy this gives sqrt(beta_i beta_j) sin(dmu).
    """

    def calculate(self, tw_init=None):
        energy = tw_init.E if tw_init is not None else 0.
        coupled, longitudinal = self._optics_coupling(energy)
        if not coupled:
            seed = tw_init
            if seed is None or (seed.beta_x == 0 and seed.beta_y == 0):
                seed = Twiss(E=energy, beta_x=1., beta_y=1.)
            self.optical_func_params(seed)
        if coupled or longitudinal:
            self.resp = self._linear_response(energy)
            return self.resp

        # In separable optics a finite Ocelot corrector is equivalent to a
        # thin angular kick at its centre.
        correctors = [*self.hcors, *self.vcors]
        nbpm = len(self.bpms)
        self.resp = np.zeros((2 * nbpm, len(correctors)))
        for j, cor in enumerate(correctors):
            for i, bpm in enumerate(self.bpms):
                if cor.lat_inx < bpm.lat_inx:
                    self.resp[i, j] = (np.cos(cor.tilt) * np.sqrt(cor.beta_x * bpm.beta_x)
                                       * np.sin(bpm.phi_x - cor.phi_x))
                    self.resp[i + nbpm, j] = (np.sin(cor.tilt) * np.sqrt(cor.beta_y * bpm.beta_y)
                                              * np.sin(bpm.phi_y - cor.phi_y))
        return self.resp


class LinacSimRM(MeasureResponseMatrix):

    def __init__(self, lattice, hcors, vcors, bpms):
        super(LinacSimRM, self).__init__(lattice, hcors, vcors, bpms)

    def calculate(self, tw_init=None):
        """Track the orbit response to each corrector at entrance reference energy tw_init.E."""
        if tw_init is None:
            raise ValueError("LinacSimRM requires tw_init with the lattice entrance energy")

        def read_bpms(write2bpms):
            return self.read_virtual_orbit(p_init=Particle(E=tw_init.E), write2bpms=write2bpms)

        return self._calculate_kick_response(read_bpms, angle_step=0.0001)


class LinacRmatrixRM(MeasureResponseMatrix):
    """First-order orbit response from corrector kicks to BPM centres in both planes."""

    def calculate(self, tw_init=None):
        energy = tw_init.E if tw_init is not None else 0.
        self.resp = self._linear_response(energy)
        return self.resp


class LinacDisperseSimRM(MeasureResponseMatrix):

    def __init__(self, lattice, hcors, vcors, bpms):
        super(LinacDisperseSimRM, self).__init__(lattice, hcors, vcors, bpms)

    def calculate(self, tw_init=None):
        """Track the change in BPM dispersion per radian of corrector kick.

        tw_init.E is the reference energy at the lattice entrance [GeV].
        Each column uses a +0.5 mrad increment and read_virtual_dispersion.
        Rows are all horizontal BPM readings followed by all vertical readings.
        Columns are horizontal correctors followed by vertical correctors.
        """
        if tw_init is None:
            raise ValueError("LinacDisperseSimRM requires tw_init with the lattice entrance energy")

        def read_bpms(write2bpms):
            return self.read_virtual_dispersion(E0=tw_init.E, write2bpms=write2bpms)

        return self._calculate_kick_response(read_bpms, angle_step=0.0005)


class LinacDisperseTmatrixRM(MeasureResponseMatrix):
    """Derivative of the BPM orbit response with respect to entrance Particle.p.

    Propagate the reference trajectory, its entrance-energy derivative and all
    corrector derivatives together. SecondTM supplies the magnetic Taylor
    tensor and RF cavities use their longitudinal derivatives. An active
    linear TransferMap has zero curvature. For other nonlinear tracking
    methods, differentiate their local map in phase space.
    No corrector settings are perturbed. Rows are X then Y, columns H then V.
    """

    def calculate(self, tw_init=None):
        if tw_init is None:
            raise ValueError("LinacDisperseTmatrixRM requires tw_init with the lattice entrance energy")
        energy = tw_init.E
        orbit = np.zeros(6)
        dispersion = np.array([0., 0., 0., 0., 0., 1.])
        kicks = np.zeros((6, len(self.hcors) + len(self.vcors)))
        mixed = np.zeros_like(kicks)
        nbpm = len(self.bpms)
        response = np.zeros((2 * nbpm, kicks.shape[1]))
        for _, maps, kick, row in self._response_sections(first_order_only=False):
            for tm in maps:
                J, dJ = map_derivatives(tm, energy, orbit, dispersion)
                mixed = J @ mixed + dJ @ kicks
                kicks = J @ kicks
                dispersion = J @ dispersion
                orbit = tm.map_function(orbit[:, None], energy)[:, 0]
                energy += tm.get_delta_e()
            if kick is not None:
                column, unit_kick = kick
                kicks[:, column] += unit_kick
            if row is not None:
                response[[row, row + nbpm]] = mixed[[0, 2]]
        self.resp = response
        return response


def _cavity_exit_reference(entrance: float, nominal_exit: float, scanned_exit: float,
                           length: float) -> tuple[float, float, float]:
    """Convert scanned RF output to the nominal momentum and reference clock."""
    energies = np.array([entrance, nominal_exit, scanned_exit])
    p0, pn, ps = np.sqrt((energies - m_e_GeV) * (energies + m_e_GeV))

    def clock_slip(exit_energy, exit_momentum):
        # For a linear energy ramp, c*t/L = (E0 + E1)/(p0 + p1).
        # Evaluate c*t/L - 1 accurately for relativistic beams (c = 1).
        numerator = m_e_GeV**2 * (1 / (entrance + p0) + 1 / (exit_energy + exit_momentum))
        return numerator / (p0 + exit_momentum)

    dtau = length * (clock_slip(scanned_exit, ps) - clock_slip(nominal_exit, pn))
    return ps / pn, (scanned_exit - nominal_exit) / pn, dtau


class LinacCavityScanRM(MeasureResponseMatrix):
    """Dispersion and its corrector response for a finite RF voltage scan.

    cavity_weights maps each selected Cavity or TWCavity to its voltage increment
    per unit station increment. voltage_step [GV] is the signed end-to-end span
    starting at the current settings. Both orbit and response differences are
    divided by voltage_step / normalization_energy [GeV], as for an electron
    beam measurement. tw_init.E supplies the independent entrance energy [GeV].
    Magnet maps and kick calibration use the nominal lattice reference throughout.
    """

    def __init__(self, lattice, hcors, vcors, bpms, *, cavity_weights, voltage_step, normalization_energy):
        super().__init__(lattice, hcors, vcors, bpms)
        if not np.isfinite(voltage_step) or voltage_step == 0.:
            raise ValueError("voltage_step must be finite and nonzero [GV]")
        if not np.isfinite(normalization_energy) or normalization_energy <= 0.:
            raise ValueError("normalization_energy must be finite and positive [GeV]")
        weights = dict(cavity_weights)
        if not weights or not all(np.isfinite(w) for w in weights.values()) or not any(weights.values()):
            raise ValueError("cavity_weights must include a nonzero weight and all weights must be finite")
        for cavity in weights:
            if not isinstance(cavity, (Cavity, TWCavity)):
                raise TypeError("Each RF scan element must be a Cavity or TWCavity")
            lattice.resolve_element_index(cavity)
        self.cavity_weights = weights
        self.voltage_step = voltage_step
        self.normalization_energy = normalization_energy

    def _endpoint(self, energy, scanned_maps, with_response):
        nbpm = len(self.bpms)
        ncor = len(self.hcors) + len(self.vcors) if with_response else 0
        readings = np.empty((4, nbpm))  # x, y, reference energy, particle delta
        response = np.empty((2 * nbpm, ncor))
        orbit = np.zeros(6)
        kicks = np.zeros((6, ncor))
        for index, nominal_maps, kick, row in self._response_sections(first_order_only=False):
            entrance = energy
            for tm in scanned_maps.get(index, nominal_maps):
                exit_energy = energy + tm.get_delta_e()
                if not np.isfinite(exit_energy) or exit_energy <= m_e_GeV:
                    raise ValueError("RF scan requires reference energies above the electron rest energy")
                if ncor:
                    J, _ = map_derivatives(tm, energy, orbit)
                    kicks = J @ kicks
                orbit = tm.map_function(orbit[:, None], energy)[:, 0]
                energy = exit_energy
            if index in scanned_maps:
                nominal = entrance + sum(tm.get_delta_e() for tm in nominal_maps)
                scale, delta, dtau = _cavity_exit_reference(entrance, nominal, energy, self.lat.sequence[index].l)
                orbit[[1, 3, 5]] *= scale
                kicks[[1, 3, 5]] *= scale
                orbit[5] += delta
                orbit[4] += dtau
                energy = nominal
            if kick is not None and ncor:
                column, unit_kick = kick
                kicks[:, column] += unit_kick
            if row is not None:
                readings[:, row] = orbit[0], orbit[2], energy, orbit[5]
                response[[row, row + nbpm]] = kicks[[0, 2]]
        return readings, response

    def _scan(self, energy, with_response, write2bpms):
        if not np.isfinite(energy) or energy <= m_e_GeV:
            raise ValueError("RF scan requires a lattice entrance energy above the electron rest energy")
        scanned_maps = {}
        for cavity, weight in self.cavity_weights.items():
            index = self.lat.resolve_element_index(cavity)
            if weight:
                atom = copy.copy(cavity.element)
                atom.v += weight * self.voltage_step
                scanned_maps[index] = [type(tm).from_element(atom, tm_type=tm.tm_type, delta_l=tm.delta_length)
                                       for tm in cavity.tms]
        baseline, nominal_response = self._endpoint(energy, {}, with_response)
        changed, scanned_response = self._endpoint(energy, scanned_maps, with_response)
        span = self.voltage_step / self.normalization_energy
        dispersion = (changed[:2] - baseline[:2]) / span
        if write2bpms:
            for i, bpm in enumerate(self.bpms):
                bpm.x, bpm.y, bpm.E, bpm.p = baseline[:, i]
                bpm.Dx, bpm.Dy = dispersion[:, i]
        return dispersion, (scanned_response - nominal_response) / span

    def read_virtual_dispersion(self, E0: float, write2bpms: bool = True) -> tuple[np.ndarray, np.ndarray]:
        """Return horizontal and vertical BPM dispersion [m] for this RF scan."""
        (dx, dy), _ = self._scan(E0, with_response=False, write2bpms=write2bpms)
        return dx, dy

    def calculate(self, tw_init=None):
        """Return the scan dispersion change per radian and store baseline BPM readings."""
        if tw_init is None:
            raise ValueError("LinacCavityScanRM requires tw_init with the lattice entrance energy")
        _, self.resp = self._scan(tw_init.E, with_response=True, write2bpms=True)
        return self.resp


class ResponseMatrixJSON:
    def __init__(self, method=None):
        self.cor_names = []
        self.bpm_names = []
        self.matrix = []
        self.method = method
        self.mode = "radian"  # or "ampere"

        self.tw_init = None  # for self.run()
        self.filename = None  # for self.run()

    def calculate(self, tw_init=None):
        """
        rewrites cor_name, bpm_name and matrix

        :param method:
        :return:
        """
        if self.method != None:
            hcors = self.method.hcors
            vcors = self.method.vcors
            bpms = self.method.bpms
            cor_names = np.append([cor.id for cor in hcors], [cor.id for cor in vcors])
            bpm_names = [bpm.id for bpm in bpms]
            matrix = self.method.calculate(tw_init=tw_init)
            self.cor_names, self.bpm_names, self.matrix = cor_names, bpm_names, matrix
        else:
            print("ResponseMatrix.method = None, Add the method, e.g. MeasureResponseMatrix")

    def get_matrix(self):
        return self.matrix

    def extract(self, cor_list, bpm_list):
        cor_list = np.array(cor_list)
        bpm_list = np.array(bpm_list)
        # print("EXTRACT given: ", cor_list)
        # print("EXTRACT given: ", bpm_list)
        cors1 = np.array(self.cor_names)
        cors2 = cor_list
        bpms1 = np.array(self.bpm_names)
        bpms2 = bpm_list
        nb1 = len(bpms1)
        # print("EXTRACT self: ", cors1)
        # print("EXTRACT self: ", bpms1)
        c_names = cors1[np.isin(cors1, cors2)]

        c_i1 = np.where(np.isin(cors1, cors2))[0]
        c_i2 = np.where(np.isin(cors2, c_names))[0]
        b_names = bpms1[np.isin(bpms1, bpms2)]
        b_i1 = np.where(np.isin(bpms1, bpms2))[0]
        b_i2 = np.where(np.isin(bpms2, b_names))[0]

        if not np.array_equal(c_names, cor_list):
            print(" Origin response matrix has no correctors:")
            print(cors2[np.isin(cors2, c_names, invert=True)])
        if not np.array_equal(b_names, bpm_list):
            print(" Origin response matrix has no BPMs:")
            print(bpm_list[b_i2[:]])

        extr_matrix = np.zeros((len(b_names) * 2, len(c_names)))
        for n in range(2):
            for i, c in enumerate(c_names):
                for j, b in enumerate(b_names):
                    x1 = self.matrix[b_i1[j] + nb1 * n, c_i1[i]]
                    extr_matrix[j + n * len(b_names), i] = x1

        return extr_matrix

    def inject(self, cor_list, bpm_list, inj_matrix):
        """
        Update some elements of the response matrix

        :param cor_list:
        :param bpm_list:
        :param inj_matrix:
        :return:
        """
        cor_list = np.array(cor_list)
        bpm_list = np.array(bpm_list)
        # print("EXTRACT given: ", cor_list)
        # print("EXTRACT given: ", bpm_list)
        cors1 = np.array(self.cor_names)
        cors2 = cor_list
        bpms1 = np.array(self.bpm_names)
        bpms2 = bpm_list
        nb1 = len(bpms1)
        # print("EXTRACT self: ", cors1)
        # print("EXTRACT self: ", bpms1)
        c_names = cors1[np.isin(cors1, cors2)]

        c_i1 = np.where(np.isin(cors1, cors2))[0]
        c_i2 = np.where(np.isin(cors2, c_names))[0]
        b_names = bpms1[np.isin(bpms1, bpms2)]
        b_i1 = np.where(np.isin(bpms1, bpms2))[0]
        b_i2 = np.where(np.isin(bpms2, b_names))[0]

        if not np.array_equal(c_names, cor_list):
            logger.warning(" ResponseMatrix.inject: Origin response matrix has no correctors:")
            print(cors2[np.isin(cors2, c_names, invert=True)])
        if not np.array_equal(b_names, bpm_list):
            logger.warning(" ResponseMatrix.inject: Origin response matrix has no BPMs:")
            print(bpm_list[b_i2[:]])

        # extr_matrix = np.zeros((len(b_names)*2, len(c_names)))

        for n in range(2):
            for i, c in enumerate(c_names):
                for j, b in enumerate(b_names):
                    # print("old matrix elem = ",b_i1[j] + nb1*n,  c_i1[i], self.matrix[b_i1[j] + nb1*n, c_i1[i]],
                    #      "inj_mat_elem =", j + n*len(b_names), i, inj_matrix[j + n*len(b_names), i])
                    self.matrix[b_i1[j] + nb1 * n, c_i1[i]] = inj_matrix[j + n * len(b_names), i]
                    # extr_matrix[j + n*len(b_names), i] = x1

        return self.matrix

    def dump(self, filename):
        dict_rmatrix = {}
        dict_rmatrix["cor_names"] = list(self.cor_names)
        dict_rmatrix["bpm_names"] = list(self.bpm_names)
        dict_rmatrix["matrix"] = list(self.matrix.flatten())

        dict_rmatrix["method_name"] = self.method.__class__.__name__ if self.method != None else "None"
        dict_rmatrix["mode"] = self.mode

        directory = os.path.dirname(filename)
        if not os.path.exists(directory):
            os.makedirs(directory)

        with open(filename, 'w+') as f:
            json.dump(dict_rmatrix, f)

    def load(self, filename):
        with open(filename, 'r') as f:
            dict_rmatrix = json.load(f)
        self.cor_names = dict_rmatrix["cor_names"]
        self.bpm_names = dict_rmatrix["bpm_names"]
        self.method_name = dict_rmatrix["method_name"]
        r_matrix = np.array(dict_rmatrix["matrix"])
        self.matrix = r_matrix.reshape(2 * len(self.bpm_names), len(self.cor_names))
        self.mode = dict_rmatrix["mode"]
        return 1

    def compare(self, rmatrix, absolut=0.001, relative=0.1):
        cors1 = np.array(self.cor_names)
        cors2 = np.array(rmatrix.cor_names)
        bpms1 = np.array(self.bpm_names)
        bpms2 = np.array(rmatrix.bpm_names)
        nb1 = len(bpms1)
        nb2 = len(bpms2)
        # c_names = np.intersect1d(cors1, cors2)
        c_names = cors1[np.isin(cors1, cors2)]
        c_i1 = np.where(np.isin(cors1, cors2))[0]
        c_i2 = np.where(np.isin(cors2, cors1))[0]
        # b_names = np.intersect1d(bpms1, bpms2)
        b_names = bpms1[np.isin(bpms1, bpms2)]
        b_i1 = np.where(np.isin(bpms1, bpms2))[0]
        b_i2 = np.where(np.isin(bpms2, bpms1))[0]
        plane = ["X", "Y"]
        for n in range(2):
            print("****************   ", plane[n], "   ****************")
            counter = 0
            for i, c in enumerate(c_names):
                for j, b in enumerate(b_names):
                    # print b_i1[j],  nb1*n, c_i1[i]
                    x1 = self.matrix[b_i1[j] + nb1 * n, c_i1[i]]
                    x2 = rmatrix.matrix[b_i2[j] + nb2 * n, c_i2[i]]
                    if abs(x1 - x2) < absolut:
                        continue
                    if abs(x1 - x2) / max(np.abs([x1, x2])) < relative:
                        continue
                    l_x1 = len(str(x1))
                    print(plane[n], c, " " * (10 - len(c)), b, " " * (10 - len(b)), "r1: ", x1, " " * (18 - l_x1),
                          "r2: ", x2)
                    counter += 1
            print("shown", counter, "elements of", len(c_names) * len(b_names))

    def show(self, list_cor=None, list_bpm=None):
        print(" " * 10, )
        for bpm in self.bpm_names:
            print(bpm, )
        print()
        for i in range(np.shape(self.matrix)[1]):
            print(self.cor_names[i] + " " * (10 - len(self.cor_names[i])), )
            # print np.array_str(self.matrix[:, i], precision=2, suppress_small=True)
            for j in range(np.shape(self.matrix)[0]):
                print("%.2f" % self.matrix[j, i], )
            print()


class ResponseMatrix:
    def __init__(self, method=None):
        self.cor_names = []
        self.bpm_names = []
        self.matrix = []
        self.method = method
        self.mode = "radian"  # or "ampere"
        self.df = None
        self.tw_init = None  # for self.run()
        self.filename = None  # for self.run()

    def bpm2x_name(self, bpm_id):
        """
        Transform bpm id to a name how it use in a control system to get horizontal beam position
        :param bpm_id:
        :return: channel for X position
        """
        return bpm_id + ".X"

    def bpm2y_name(self, bpm_id):
        """
        Transform bpm id to a name how it use in a control system to get vertical beam position

        :param bpm_id:
        :return: channel for Y position
        """
        return bpm_id + ".Y"

    def xy_names2bpm_id(self, xy_names):
        """
        transform BPM channels to bpm ids

        :param xy_names:
        :return:
        """
        bpm_ids = [bpm.replace(".X", "") for bpm in xy_names if ".X" in bpm]
        return bpm_ids

    def calculate(self, tw_init=None):
        """
        rewrites cor_name, bpm_name and matrix

        :param method:
        :return:
        """
        if self.method != None:
            hcors = self.method.hcors
            vcors = self.method.vcors
            bpms = self.method.bpms
            cor_names = np.append([cor.id for cor in hcors], [cor.id for cor in vcors])
            bpm_names = [bpm.id for bpm in bpms]
            matrix = self.method.calculate(tw_init=tw_init)
            self.df = self.data2df(matrix=matrix, bpm_names=bpm_names, cor_names=cor_names)
            self.cor_names, self.bpm_names, self.matrix = cor_names, bpm_names, matrix
        else:
            print("ResponseMatrix.method = None, Add the method, e.g. MeasureResponseMatrix")

    def get_matrix(self):
        return self.matrix

    def extract_df_slice(self, cor_list, bpm_list):
        bpm_x = [self.bpm2x_name(bpm) for bpm in bpm_list]
        bpm_y = [self.bpm2y_name(bpm) for bpm in bpm_list]
        rows = bpm_x + bpm_y
        cols = list(cor_list)
        cor_list_exist = np.array([item in self.df.columns for item in cor_list])
        bpm_list_exist = np.array([item in self.df.index for item in rows])

        if all(cor_list_exist) and all(bpm_list_exist):
            df_slice = self.df.loc[rows, cols]
            return df_slice
        else:
            print("correctors are not in the RM")
            print(np.array(cor_list)[~cor_list_exist])
            print()
            print("BPMs are not in the RM")
            print(np.array(rows)[~bpm_list_exist])
            return None

    def extract(self, cor_list, bpm_list):
        df_slice = self.extract_df_slice(cor_list, bpm_list)
        # pandas 3 exposes shared arrays as read-only; callers edit this slice
        # before passing it back to inject().
        return df_slice.to_numpy(copy=True)

    def retrieve_from_scan(self, df_scan):
        from sklearn.linear_model import LinearRegression

        bpm_x = [self.bpm2x_name(bpm) for bpm in self.bpm_names]
        bpm_y = [self.bpm2y_name(bpm) for bpm in self.bpm_names]
        bpm_names_xy = bpm_x + bpm_y

        x = df_scan.loc[:, self.cor_names].values
        y = df_scan.loc[:, bpm_names_xy].values

        reg = LinearRegression().fit(x, y)
        self.df = self.data2df(matrix=reg.coef_, bpm_names=self.bpm_names, cor_names=self.cor_names)
        self.matrix = self.df.to_numpy(copy=True)
        return self.df

    def clean_rm(self, coupling=True):
        if self.method != None:
            ring = isinstance(self.method, RingRM)
            for correctors, cross_plane in ((self.method.hcors, self.bpm2y_name),
                                             (self.method.vcors, self.bpm2x_name)):
                for cor in correctors:
                    for bpm in self.method.bpms:
                        if not ring and bpm.s < cor.s:
                            self.df.loc[[self.bpm2x_name(bpm.id), self.bpm2y_name(bpm.id)], cor.id] = 0
                        elif not coupling:
                            self.df.loc[cross_plane(bpm.id), cor.id] = 0
            self.matrix = self.df.to_numpy(copy=True)

        else:
            print("ResponseMatrix.method = None, Add the method, e.g. MeasureResponseMatrix")

    def inject(self, cor_list, bpm_list, inj_matrix):
        """
        Update some elements of the response matrix

        :param cor_list:
        :param bpm_list:
        :param inj_matrix:
        :return:
        """
        df_slice = self.data2df(matrix=inj_matrix, bpm_names=bpm_list, cor_names=cor_list)
        self.df.update(df_slice)
        self.matrix = self.df.to_numpy(copy=True)
        return self.matrix

    def data2df(self, matrix, bpm_names, cor_names):
        bpm_x = [self.bpm2x_name(bpm) for bpm in bpm_names]
        bpm_y = [self.bpm2y_name(bpm) for bpm in bpm_names]
        df = pd.DataFrame(matrix, columns=cor_names, index=bpm_x + bpm_y)
        return df

    def df2data(self):
        self.cor_names = list(self.df.columns.values)
        bpms_all = list(self.df.index.values)
        self.bpm_names = self.xy_names2bpm_id(bpms_all)
        self.matrix = self.df.to_numpy(copy=True)

    def dump(self, filename):
        df = self.data2df(matrix=self.matrix, bpm_names=self.bpm_names, cor_names=self.cor_names)
        directory = os.path.dirname(filename)
        if directory != "" and not os.path.exists(directory):
            os.makedirs(directory)
        df.to_pickle(filename)

    def load(self, filename):
        self.df = pd.read_pickle(filename)
        self.df2data()
        return 1

    def compare(self, rmatrix, absolut=0.001, relative=0.1):
        cors1 = np.array(self.cor_names)
        cors2 = np.array(rmatrix.cor_names)
        bpms1 = np.array(self.bpm_names)
        bpms2 = np.array(rmatrix.bpm_names)
        nb1 = len(bpms1)
        nb2 = len(bpms2)
        # c_names = np.intersect1d(cors1, cors2)
        c_names = cors1[np.isin(cors1, cors2)]
        c_i1 = np.where(np.isin(cors1, cors2))[0]
        c_i2 = np.where(np.isin(cors2, cors1))[0]
        # b_names = np.intersect1d(bpms1, bpms2)
        b_names = bpms1[np.isin(bpms1, bpms2)]
        b_i1 = np.where(np.isin(bpms1, bpms2))[0]
        b_i2 = np.where(np.isin(bpms2, bpms1))[0]
        plane = ["X", "Y"]
        for n in range(2):
            print("****************   ", plane[n], "   ****************")
            counter = 0
            for i, c in enumerate(c_names):
                for j, b in enumerate(b_names):
                    # print b_i1[j],  nb1*n, c_i1[i]
                    x1 = self.matrix[b_i1[j] + nb1 * n, c_i1[i]]
                    x2 = rmatrix.matrix[b_i2[j] + nb2 * n, c_i2[i]]
                    if abs(x1 - x2) < absolut:
                        continue
                    if abs(x1 - x2) / max(np.abs([x1, x2])) < relative:
                        continue
                    l_x1 = len(str(x1))
                    print(plane[n], c, " " * (10 - len(c)), b, " " * (10 - len(b)), "r1: ", x1, " " * (18 - l_x1),
                          "r2: ", x2)
                    counter += 1
            print("shown", counter, "elements of", len(c_names) * len(b_names))

    def show(self, list_cor=None, list_bpm=None):
        print(" " * 10, )
        for bpm in self.bpm_names:
            print(bpm, )
        print()
        for i in range(np.shape(self.matrix)[1]):
            print(self.cor_names[i] + " " * (10 - len(self.cor_names[i])), )
            # print np.array_str(self.matrix[:, i], precision=2, suppress_small=True)
            for j in range(np.shape(self.matrix)[0]):
                print("%.2f" % self.matrix[j, i], )
            print()
