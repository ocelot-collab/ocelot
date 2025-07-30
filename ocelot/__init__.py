"""
general ocelot description
"""

__version__ = '25.07.1'


__all__ = ['Twiss', "Beam", "Particle", "get_current", "get_envelope", "generate_parray",           # beam
            "ellipse_from_twiss", "ParticleArray",  "global_slice_analysis", 'gauss_from_twiss',    # beam

            "save_particle_array", "load_particle_array",                                           # io

            'fodo_parameters', 'lattice_transfer_map', "Navigator", 'twiss', "MethodTM",            # optics

            'Element', 'Multipole', 'Quadrupole', 'RBend', "Matrix", "UnknownElement",              # elements
            'SBend', 'Bend', 'Drift', 'Undulator', 'Hcor', "Solenoid", "TDCavity",                  # elements
            'Vcor', "Sextupole", "Monitor", "Marker", "Octupole", "Cavity",  "Aperture",            # elements

            "match", "match_tunes",                                                                 # match

            "tracking_step", "create_track_list", "track_nturns", "freq_analysis",                  # track
            "contour_da", "track_nturns_mpi", "nearest_particle", "stable_particles",               # track
            "spectrum", "track",                                                                    # track
            "pi", "m_e_eV", "m_e_MeV", "m_e_GeV", "speed_of_light",                             # globals
            "compensate_chromaticity",                                                          # chromaticity
            "EbeamParams",                                                                      # beam_params
            "CSR",                                                                              # csr
            "SpaceCharge", "LSC",                                                               # sc
            "Wake", "WakeTable", "WakeKick", "WakeTableDechirperOffAxis",                       # wake
            "BeamTransform", "SmoothBeam", "EmptyProc", "PhysProc", "LaserHeater",
            "LaserModulator", "SpontanRadEffects", "PhaseSpaceAperture",
            "RectAperture", "EllipticalAperture",
            "CopyBeam", "SaveBeam", "LatticeEnergyProfile",
            "MagneticLattice", "merger",                                                        # magnetic_lattice
            "np",                                                                               # numpy

           "CavityTM", "TransferMap",                                                   # transformations
           "KickTM", "MultipoleTM", "PulseTM", "RungeKuttaTM", "RungeKuttaTrTM",        # transformations
           "SecondTM", "TWCavityTM", "UndulatorTestTM", "TMTypes",                      # transformations

            "CavityParams", "FirstOrderParams" ,"KickParams", "MultipoleParams",        # tm parameters
            "RungeKuttaParams", "SecondOrderParams", "UndulatorTestParams"              # tm parameters

           ]

import numpy as np
from ocelot.cpbd.magnetic_lattice import MagneticLattice, merger

from ocelot.cpbd.beam import *
from ocelot.cpbd.tm_params import *
from ocelot.cpbd.transformations import *
from ocelot.cpbd.elements import *
from ocelot.cpbd.match import *
from ocelot.cpbd.track import *
from ocelot.common.globals import *
from ocelot.common.ocelog import *
from ocelot.cpbd.chromaticity import *
from ocelot.cpbd.beam_params import *
from ocelot.cpbd.io import *
from ocelot.cpbd.sc import *
from ocelot.cpbd.csr import *
from ocelot.cpbd.wake3D import *
from ocelot.cpbd.physics_proc import *
from ocelot.cpbd.navi import Navigator


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

