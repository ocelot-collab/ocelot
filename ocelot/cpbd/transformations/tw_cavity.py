from ocelot.cpbd.transformations.cavity import CavityTM


class TWCavityTM(CavityTM):
    """Track with travelling-wave transverse optics and the shared RF body kick.

    TWCavityAtom supplies the edge lenses and adiabatic transverse transport.
    The longitudinal tracking algorithm and slice-voltage convention are shared
    with the standing-wave cavity.
    """
