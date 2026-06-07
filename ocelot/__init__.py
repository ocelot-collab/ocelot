"""
general ocelot description
"""

__version__ = '26.06.0'


__all__ = [
    # === Beam & Particle Physics ===
    'Twiss', "Beam", "Particle", "get_current", "get_envelope", "generate_parray",
    "ellipse_from_twiss", "ParticleArray", "global_slice_analysis", 'gauss_from_twiss',

    # === Input/Output ===
    "save_particle_array", "load_particle_array",

    # === Optics & Navigation ===
    'fodo_parameters', 'lattice_transfer_map', "Navigator", 'twiss', "MethodTM",

    # === Lattice Elements ===
    'Element', 'Multipole', 'Quadrupole', 'RBend', "Matrix", "UnknownElement",
    'SBend', 'Bend', 'Drift', 'Undulator', 'Hcor', "Solenoid", "TDCavity",
    'Vcor', "Sextupole", "Monitor", "Marker", "Octupole", "Cavity", "Aperture",

    # === Lattice Matching ===
    "match", "match_tunes",

    # === Tracking ===
    "tracking_step", "create_track_list", "track_nturns", "freq_analysis",
    "contour_da", "track_nturns_mpi", "nearest_particle", "stable_particles",
    "spectrum", "track",

    # === Global Constants ===
    "pi", "m_e_eV", "m_e_MeV", "m_e_GeV", "speed_of_light",

    # === Beam Dynamics & Phenomena ===
    "compensate_chromaticity", "EbeamParams", "CSR", "SpaceCharge", "LSC",

    # === Wake Effects & Physics Processes ===
    "Wake", "WakeTable", "WakeKick", "WakeTableDechirperOffAxis", "LongWake", "LinLongWake",
    "BeamTransform", "SmoothBeam", "EmptyProc", "PhysProc", "LaserHeater",
    "LaserModulator", "SpontanRadEffects", "PhaseSpaceAperture",
    "RectAperture", "EllipticalAperture", "CopyBeam", "SaveBeam", "LatticeEnergyProfile",

    # === Magnetic Lattice ===
    "MagneticLattice", "merger",

    # === External Dependencies ===
    "np",

    # === Transfer Maps & Transformations ===
    "CavityTM", "TransferMap", "KickTM", "MultipoleTM", "PulseTM",
    "RungeKuttaGlobalTM", "RungeKuttaOcelotTM", "RungeKuttaTM", "RungeKuttaTrTM", "SecondTM", "TWCavityTM",
    "UndulatorTestTM", "TMTypes",

    # === Transfer Map Parameters ===
    "CavityParams", "FirstOrderParams", "KickParams", "MultipoleParams",
    "RungeKuttaParams", "SecondOrderParams", "UndulatorTestParams"
]

# ============================================================================
# External Dependencies
# ============================================================================
import numpy as np

# ============================================================================
# Lattice Infrastructure
# ============================================================================
from ocelot.cpbd.magnetic_lattice import MagneticLattice, merger
from ocelot.cpbd.navi import Navigator

# ============================================================================
# Beam & Particle Physics
# ============================================================================
from ocelot.cpbd.beam.core import Twiss, Beam
from ocelot.cpbd.beam.particle import ParticleArray, Particle
from ocelot.cpbd.beam.beam import ellipse_from_twiss, gauss_from_twiss
from ocelot.cpbd.beam.generator import generate_parray
from ocelot.cpbd.beam.analysis import get_current, get_envelope, global_slice_analysis

# ============================================================================
# Input/Output
# ============================================================================
from ocelot.cpbd.io import *

# ============================================================================
# Optics & Matching
# ============================================================================
from ocelot.cpbd.match import *

# ============================================================================
# Transfer Maps & Parameters
# ============================================================================
from ocelot.cpbd.tm_params import *
from ocelot.cpbd.transformations import *

# ============================================================================
# Lattice Elements
# ============================================================================
from ocelot.cpbd.elements import *

# ============================================================================
# Tracking & Analysis
# ============================================================================
from ocelot.cpbd.track import *

# ============================================================================
# Global Constants & Utilities
# ============================================================================
from ocelot.common.globals import *
from ocelot.common.ocelog import *

# ============================================================================
# Beam Dynamics & Phenomena
# ============================================================================
from ocelot.cpbd.chromaticity import *
from ocelot.cpbd.beam_params import *
from ocelot.cpbd.csr import *
from ocelot.cpbd.sc import *
from ocelot.cpbd.wake3D import *
from ocelot.cpbd.physics_proc import *


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
