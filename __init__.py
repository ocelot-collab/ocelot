'''
general ocelot description
'''

__version__ = '16.8.rc'

__all__ = ['Twiss', 'twiss', "Beam", "Particle", "get_current", "get_envelope",  # beam
            "ellipse_from_twiss", "ParticleArray",  # beam
           'fodo_parameters', 'lattice_transfer_map', 'TransferMap', 'gauss_from_twiss',  # optics
           "get_map", "MethodTM", "SecondTM", "KickTM", "CavityTM", "UndulatorTestTM",  # optics
           'Element', 'Multipole', 'Quadrupole', 'RBend', "Matrix", "UnknownElement",  # elements
           'SBend', 'Bend', 'Drift', 'Undulator', 'Hcor',  # elements
           'Vcor', "Sextupole", "Monitor", "Marker", "Octupole", "Cavity", "Edge",  # elements
           "Sequence", "Solenoid",  # elements
           "match", "match_tunes",  # match
           "Navigator", "tracking_step", "create_track_list", "track_nturns", "freq_analysis",  # track
            "contour_da", "track_nturns_mpi", "nearest_particle", "stable_particles",  # track
            "spectrum", "track",  # track
           "pi", "m_e_eV", "m_e_MeV", "m_e_GeV",  # globals
           "compensate_chromaticity",  # chromaticity
           "EbeamParams",  # e_beam_params
           "write_lattice",  # io
           "CSR", "SpaceCharge", "Wake", "WakeTable",
           "MagneticLattice",
           ]

from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.beam import *
from ocelot.cpbd.optics import *
from ocelot.cpbd.elements import *
from ocelot.cpbd.match import *
from ocelot.cpbd.track import *
from ocelot.common.globals import *
from ocelot.common.logging import Logger
from ocelot.cpbd.chromaticity import *
from ocelot.cpbd.e_beam_params import *
from ocelot.cpbd.io import *
from ocelot.cpbd.sc import *
from ocelot.cpbd.csr import *
from ocelot.cpbd.wake3D import *
print('initializing ocelot...')
logger = Logger()
xrange=range
