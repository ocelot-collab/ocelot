"""
general ocelot description
"""

__version__ = '20.04.0'


__all__ = ['Twiss', "Beam", "Particle", "get_current", "get_envelope", "generate_parray",           # beam
            "ellipse_from_twiss", "ParticleArray",  "global_slice_analysis", 'gauss_from_twiss',    # beam

            "save_particle_array", "load_particle_array", "write_lattice",                          # io

            'fodo_parameters', 'lattice_transfer_map', 'TransferMap', "Navigator", 'twiss',    # optics
            "get_map", "MethodTM", "SecondTM", "KickTM", "CavityTM", "UndulatorTestTM",        # optics

            'Element', 'Multipole', 'Quadrupole', 'RBend', "Matrix", "UnknownElement",              # elements
            'SBend', 'Bend', 'Drift', 'Undulator', 'Hcor',  "Sequence", "Solenoid", "TDCavity",     # elements
            'Vcor', "Sextupole", "Monitor", "Marker", "Octupole", "Cavity", "Edge",  "Aperture",    # elements

            "match", "match_tunes",                                                          # match

            "tracking_step", "create_track_list", "track_nturns", "freq_analysis",           # track
            "contour_da", "track_nturns_mpi", "nearest_particle", "stable_particles",        # track
            "spectrum", "track",                                                             # track
            "pi", "m_e_eV", "m_e_MeV", "m_e_GeV", "speed_of_light",                             # globals
            "compensate_chromaticity",                                                          # chromaticity
            "EbeamParams",                                                                      # beam_params
            "CSR",                                                                              # csr
            "SpaceCharge", "LSC",                                                               # sc
            "Wake", "WakeTable", "WakeKick", "WakeTableDechirperOffAxis",                       # wake
            "BeamTransform", "SmoothBeam", "EmptyProc", "PhysProc", "LaserHeater",
            "LaserModulator", "SpontanRadEffects", "PhaseSpaceAperture",
            "MagneticLattice", "merger",            # magnetic_lattice
            "np", # numpy


           ]

import numpy as np
from ocelot.cpbd.magnetic_lattice import MagneticLattice, merger

from ocelot.cpbd.beam import ParticleArray, global_slice_analysis,Particle, Beam, Twiss, get_current, \
    get_envelope, generate_parray, ellipse_from_twiss

from ocelot.cpbd.optics import lattice_transfer_map, TransferMap, Navigator, twiss, get_map, MethodTM, \
    SecondTM, KickTM, CavityTM, UndulatorTestTM

from ocelot.cpbd.elements import *
from ocelot.cpbd.match import match, match_tunes
from ocelot.cpbd.track import *
from ocelot.common.globals import pi, m_e_eV, m_e_MeV, m_e_GeV, speed_of_light
from ocelot.common.ocelog import *
from ocelot.cpbd.chromaticity import compensate_chromaticity
from ocelot.cpbd.beam_params import EbeamParams
from ocelot.cpbd.io import save_particle_array, load_particle_array, write_lattice
from ocelot.cpbd.sc import SpaceCharge, LSC
from ocelot.cpbd.csr import CSR
from ocelot.cpbd.wake3D import Wake, WakeTable, WakeKick, WakeTableDechirperOffAxis

from ocelot.cpbd.physics_proc import BeamTransform, SmoothBeam, EmptyProc, PhysProc, LaserHeater, \
    LaserModulator, SpontanRadEffects, PhaseSpaceAperture


print('initializing ocelot...')

try:
    import numba
except:
    print("import: module NUMBA is not installed. Install it to speed up calculation")

try:
    import pyfftw
except:
    print("import: module PYFFTW is not installed. Install it to speed up calculation")

try:
    import numexpr
except:
    print("import: module NUMEXPR is not installed. Install it to speed up calculation")

