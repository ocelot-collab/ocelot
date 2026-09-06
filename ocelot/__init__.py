"""
Public Ocelot package facade.

The root package intentionally stays lightweight. Public names are loaded on
first access so that ``import ocelot`` and narrow submodule imports do not pay
for the full accelerator, radiation, I/O, and plotting stack.
"""

from importlib import import_module

__version__ = '26.06.1'


class _CopyCompat:
    """Compatibility object for legacy ``from ocelot import *`` workflows.

    Older examples used ``copy(obj)`` from the Ocelot root namespace, while
    some test/config modules import Python's ``copy`` module before star
    importing Ocelot and then call ``copy.copy`` or ``copy.deepcopy``.
    """

    def __call__(self, obj):
        from copy import copy

        return copy(obj)

    @staticmethod
    def copy(obj):
        from copy import copy

        return copy(obj)

    @staticmethod
    def deepcopy(obj, memo=None):
        from copy import deepcopy

        if memo is None:
            return deepcopy(obj)
        return deepcopy(obj, memo)


_COPY_COMPAT = _CopyCompat()


__all__ = [
    # === Beam & Particle Physics ===
    'Twiss', "Beam", "Particle", "get_current", "get_envelope", "generate_parray",
    "ellipse_from_twiss", "ParticleArray", "global_slice_analysis", 'gauss_from_twiss',

    # === Input/Output ===
    "save_particle_array", "load_particle_array",

    # === Optics & Navigation ===
    'fodo_parameters', 'lattice_transfer_map', "Navigator", 'twiss',
    'periodic_twiss', 'UnstableLatticeError', "MethodTM",

    # === Lattice Elements ===
    'Element', 'Multipole', 'Quadrupole', 'RBend', "Matrix", "UnknownElement",
    'SBend', 'Bend', 'Drift', 'Undulator', 'Hcor', "Solenoid", "TDCavity",
    'Vcor', "Sextupole", "Monitor", "Marker", "Octupole", "Cavity", "Aperture",
    "XYQuadrupole",

    # === Lattice Matching ===
    "match", "match_tunes",

    # === Tracking ===
    "tracking_step", "create_track_list", "track_nturns", "freq_analysis",
    "contour_da", "track_nturns_mpi", "nearest_particle", "stable_particles",
    "spectrum", "track", "lattice_track",

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
    "np", "copy", "deepcopy",

    # === Transfer Maps & Transformations ===
    "CavityTM", "TransferMap", "ExactDriftTM", "KickTM", "MultipoleTM", "PulseTM",
    "RungeKuttaGlobalTM", "RungeKuttaOcelotTM", "RungeKuttaTM", "RungeKuttaTrTM", "SecondTM", "TWCavityTM",
    "UndulatorTestTM", "TMTypes",

    # === Transfer Map Parameters ===
    "CavityParams", "ExactDriftParams", "FirstOrderParams", "KickParams", "MultipoleParams",
    "RungeKuttaParams", "SecondOrderParams", "UndulatorTestParams"
]


_LAZY_EXPORTS = {
    # External dependencies and logging compatibility.
    "np": ("numpy", None),
    "copy": (None, None),
    "deepcopy": ("copy", "deepcopy"),
    "logging": ("logging", None),
    "ocelog": ("ocelot.common.ocelog", "ocelog"),

    # Lattice infrastructure.
    "MagneticLattice": ("ocelot.cpbd.magnetic_lattice", "MagneticLattice"),
    "merger": ("ocelot.cpbd.magnetic_lattice", "merger"),
    "Navigator": ("ocelot.cpbd.navi", "Navigator"),

    # Beam and particle physics.
    "Twiss": ("ocelot.cpbd.beam.core", "Twiss"),
    "Beam": ("ocelot.cpbd.beam.core", "Beam"),
    "ParticleArray": ("ocelot.cpbd.beam.particle", "ParticleArray"),
    "Particle": ("ocelot.cpbd.beam.particle", "Particle"),
    "ellipse_from_twiss": ("ocelot.cpbd.beam.beam", "ellipse_from_twiss"),
    "gauss_from_twiss": ("ocelot.cpbd.beam.beam", "gauss_from_twiss"),
    "generate_parray": ("ocelot.cpbd.beam.generator", "generate_parray"),
    "get_current": ("ocelot.cpbd.beam.analysis", "get_current"),
    "get_envelope": ("ocelot.cpbd.beam.analysis", "get_envelope"),
    "global_slice_analysis": ("ocelot.cpbd.beam.analysis", "global_slice_analysis"),

    # Input/output.
    "save_particle_array": ("ocelot.cpbd.io", "save_particle_array"),
    "load_particle_array": ("ocelot.cpbd.io", "load_particle_array"),

    # Optics and matching.
    "fodo_parameters": ("ocelot.cpbd.optics", "fodo_parameters"),
    "lattice_transfer_map": ("ocelot.cpbd.optics", "lattice_transfer_map"),
    "twiss": ("ocelot.cpbd.optics", "twiss"),
    "periodic_twiss": ("ocelot.cpbd.optics", "periodic_twiss"),
    "UnstableLatticeError": ("ocelot.cpbd.optics", "UnstableLatticeError"),
    "MethodTM": ("ocelot.cpbd.optics", "MethodTM"),
    "match": ("ocelot.cpbd.match", "match"),
    "match_tunes": ("ocelot.cpbd.match", "match_tunes"),

    # Lattice elements.
    "Element": ("ocelot.cpbd.elements.element", "Element"),
    "Multipole": ("ocelot.cpbd.elements.multipole", "Multipole"),
    "Quadrupole": ("ocelot.cpbd.elements.quadrupole", "Quadrupole"),
    "RBend": ("ocelot.cpbd.elements.rbend", "RBend"),
    "Matrix": ("ocelot.cpbd.elements.matrix", "Matrix"),
    "UnknownElement": ("ocelot.cpbd.elements.unknown_element", "UnknownElement"),
    "SBend": ("ocelot.cpbd.elements.sbend", "SBend"),
    "Bend": ("ocelot.cpbd.elements.bend", "Bend"),
    "Drift": ("ocelot.cpbd.elements.drift", "Drift"),
    "Undulator": ("ocelot.cpbd.elements.undulator", "Undulator"),
    "Hcor": ("ocelot.cpbd.elements.hcor", "Hcor"),
    "Solenoid": ("ocelot.cpbd.elements.solenoid", "Solenoid"),
    "TDCavity": ("ocelot.cpbd.elements.tdcavity", "TDCavity"),
    "Vcor": ("ocelot.cpbd.elements.vcor", "Vcor"),
    "Sextupole": ("ocelot.cpbd.elements.sextupole", "Sextupole"),
    "Monitor": ("ocelot.cpbd.elements.monitor", "Monitor"),
    "Marker": ("ocelot.cpbd.elements.marker", "Marker"),
    "Octupole": ("ocelot.cpbd.elements.octupole", "Octupole"),
    "Cavity": ("ocelot.cpbd.elements.cavity", "Cavity"),
    "Aperture": ("ocelot.cpbd.elements.aperture", "Aperture"),
    "XYQuadrupole": ("ocelot.cpbd.elements.xyquadruple", "XYQuadrupole"),

    # Tracking and analysis.
    "tracking_step": ("ocelot.cpbd.track", "tracking_step"),
    "create_track_list": ("ocelot.cpbd.track", "create_track_list"),
    "track_nturns": ("ocelot.cpbd.track", "track_nturns"),
    "freq_analysis": ("ocelot.cpbd.track", "freq_analysis"),
    "contour_da": ("ocelot.cpbd.track", "contour_da"),
    "track_nturns_mpi": ("ocelot.cpbd.track", "track_nturns_mpi"),
    "nearest_particle": ("ocelot.cpbd.track", "nearest_particle"),
    "stable_particles": ("ocelot.cpbd.track", "stable_particles"),
    "spectrum": ("ocelot.cpbd.track", "spectrum"),
    "track": ("ocelot.cpbd.track", "track"),
    "lattice_track": ("ocelot.cpbd.track", "lattice_track"),

    # Global constants.
    "pi": ("ocelot.common.globals", "pi"),
    "m_e_eV": ("ocelot.common.globals", "m_e_eV"),
    "m_e_MeV": ("ocelot.common.globals", "m_e_MeV"),
    "m_e_GeV": ("ocelot.common.globals", "m_e_GeV"),
    "speed_of_light": ("ocelot.common.globals", "speed_of_light"),

    # Beam dynamics and phenomena.
    "compensate_chromaticity": ("ocelot.cpbd.chromaticity", "compensate_chromaticity"),
    "EbeamParams": ("ocelot.cpbd.beam_params", "EbeamParams"),
    "CSR": ("ocelot.cpbd.csr", "CSR"),
    "SpaceCharge": ("ocelot.cpbd.sc", "SpaceCharge"),
    "LSC": ("ocelot.cpbd.sc", "LSC"),

    # Wake effects and physics processes.
    "Wake": ("ocelot.cpbd.wake3D", "Wake"),
    "WakeTable": ("ocelot.cpbd.wake3D", "WakeTable"),
    "WakeKick": ("ocelot.cpbd.wake3D", "WakeKick"),
    "WakeTableDechirperOffAxis": ("ocelot.cpbd.wake3D", "WakeTableDechirperOffAxis"),
    "LongWake": ("ocelot.cpbd.wake3D", "LongWake"),
    "LinLongWake": ("ocelot.cpbd.wake3D", "LinLongWake"),
    "BeamTransform": ("ocelot.cpbd.physics_proc", "BeamTransform"),
    "SmoothBeam": ("ocelot.cpbd.physics_proc", "SmoothBeam"),
    "EmptyProc": ("ocelot.cpbd.physics_proc", "EmptyProc"),
    "PhysProc": ("ocelot.cpbd.physics_proc", "PhysProc"),
    "LaserHeater": ("ocelot.cpbd.physics_proc", "LaserHeater"),
    "LaserModulator": ("ocelot.cpbd.physics_proc", "LaserModulator"),
    "SpontanRadEffects": ("ocelot.cpbd.physics_proc", "SpontanRadEffects"),
    "PhaseSpaceAperture": ("ocelot.cpbd.physics_proc", "PhaseSpaceAperture"),
    "RectAperture": ("ocelot.cpbd.physics_proc", "RectAperture"),
    "EllipticalAperture": ("ocelot.cpbd.physics_proc", "EllipticalAperture"),
    "CopyBeam": ("ocelot.cpbd.physics_proc", "CopyBeam"),
    "SaveBeam": ("ocelot.cpbd.physics_proc", "SaveBeam"),
    "LatticeEnergyProfile": ("ocelot.cpbd.physics_proc", "LatticeEnergyProfile"),

    # Transfer maps and transformations.
    "CavityTM": ("ocelot.cpbd.transformations.cavity", "CavityTM"),
    "TransferMap": ("ocelot.cpbd.transformations.transfer_map", "TransferMap"),
    "ExactDriftTM": ("ocelot.cpbd.transformations.exact_drift", "ExactDriftTM"),
    "KickTM": ("ocelot.cpbd.transformations.kick", "KickTM"),
    "MultipoleTM": ("ocelot.cpbd.transformations.multipole", "MultipoleTM"),
    "PulseTM": ("ocelot.cpbd.transformations.pulse", "PulseTM"),
    "RungeKuttaGlobalTM": ("ocelot.cpbd.transformations.runge_kutta", "RungeKuttaGlobalTM"),
    "RungeKuttaOcelotTM": ("ocelot.cpbd.transformations.runge_kutta", "RungeKuttaOcelotTM"),
    "RungeKuttaTM": ("ocelot.cpbd.transformations.runge_kutta", "RungeKuttaTM"),
    "RungeKuttaTrTM": ("ocelot.cpbd.transformations.runge_kutta_tr", "RungeKuttaTrTM"),
    "SecondTM": ("ocelot.cpbd.transformations.second_order", "SecondTM"),
    "TWCavityTM": ("ocelot.cpbd.transformations.tw_cavity", "TWCavityTM"),
    "UndulatorTestTM": ("ocelot.cpbd.transformations.undulator_test", "UndulatorTestTM"),
    "TMTypes": ("ocelot.cpbd.transformations.transformation", "TMTypes"),

    # Transfer map parameters.
    "CavityParams": ("ocelot.cpbd.tm_params.cavity_params", "CavityParams"),
    "ExactDriftParams": ("ocelot.cpbd.tm_params.exact_drift_params", "ExactDriftParams"),
    "FirstOrderParams": ("ocelot.cpbd.tm_params.first_order_params", "FirstOrderParams"),
    "KickParams": ("ocelot.cpbd.tm_params.kick_params", "KickParams"),
    "MultipoleParams": ("ocelot.cpbd.tm_params.multipole_params", "MultipoleParams"),
    "RungeKuttaParams": ("ocelot.cpbd.tm_params.runge_kutta_params", "RungeKuttaParams"),
    "SecondOrderParams": ("ocelot.cpbd.tm_params.second_order_params", "SecondOrderParams"),
    "UndulatorTestParams": ("ocelot.cpbd.tm_params.undulator_test_params", "UndulatorTestParams"),
}


def __getattr__(name):
    try:
        module_name, attr_name = _LAZY_EXPORTS[name]
    except KeyError as exc:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc

    if name == "copy":
        value = _COPY_COMPAT
    else:
        module = import_module(module_name)
        value = module if attr_name is None else getattr(module, attr_name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__) | set(_LAZY_EXPORTS))
