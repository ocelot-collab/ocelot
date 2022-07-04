import ocelot.cpbd.elements as elements
from ocelot.cpbd.beam import Twiss
import logging

try:
    import tfs
except ImportError:
    logging.getLogger(__file__).warning(
        f"Optional package: TFS-Pandas missing.  {__name__} will lack functionality."
    )




def tfs_to_cell_and_optics(tfs_path):
    tfs_table = tfs.read(tfs_path)

    cell = convert_tfs_lattice(tfs_table)
    optics = optics_from_tfs(tfs_table)
    return cell, optics

def optics_from_tfs(tfs_table):
    twiss = Twiss()

    try:
        header = tfs_table.headers
    except AttributeError:
        tfs_table = tfs.read(tfs_table)
        header = tfs_table.headers

    twiss.E = header["ENERGY"]
    twiss.emit_x = header["EX"]
    twiss.emit_y = header["EY"]
    rel_gamma = header["GAMMA"]
    twiss.emit_xn = twiss.emit_x * rel_gamma
    twiss.emit_yn = twiss.emit_y * rel_gamma

    first = tfs_table.iloc[0]
    twiss.alpha_x = first["ALFX"]
    twiss.alpha_y = first["ALFY"]
    twiss.beta_x = first["BETX"]
    twiss.beta_y = first["BETY"]
    twiss.Dx = first["DX"]
    twiss.Dy = first["DY"]
    twiss.Dxp = first["DPX"]
    twiss.Dyp = first["DPY"]
    twiss.x = first["X"]
    twiss.y = first["Y"]
    twiss.xp = first["PX"]
    twiss.yp = first["PY"]

    return twiss


class UnsupportedMADXElementType(RuntimeError): pass

def convert_tfs_lattice(tfs_table, converter=None):
    if converter is None:
        converter = MADXLatticeConverter()

    for row in tfs_table.itertuples():
        yield converter.dispatch(row)


class MADXLatticeConverter:
    def dispatch(self, row):
        etype = row.KEYWORD
        etypelow = etype.lower()
        try:
            return getattr(self, f"make_{etypelow}")(row)
        except AttributeError:
            raise UnsupportedMADXElementType(
                f"Unsupported element type: {etype}."
            )

    def make_marker(self, row):
        return elements.Marker(eid=row.NAME)

    def make_monitor(self, row):
        return elements.Monitor(eid=row.NAME, l = row.L)

    def make_drift(self, row):
        return elements.Drift(eid=row.NAME, l=row.L)

    def make_sbend(self, row):
        return elements.SBend(eid=row.NAME,
                              l=row.L,
                              angle=row.ANGLE,
                              tilt=row.TILT,
                              e1=row.E1,
                              e2=row.E2)

    def make_rbend(self, row):
        return elements.RBend(eid=row.NAME,
                              l=row.L,
                              angle=row.ANGLE,
                              k1=row.K1L/row.L,
                              k2=row.K2L/row.L,
                              tilt=row.TILT,
                              e1=row.E1,
                              e2=row.E2,
                              fint=row.FINT,
                              fintx=row.FINTX)

    def make_quadrupole(self, row):
        return elements.Quadrupole(eid=row.NAME,
                                   l=row.L,
                                   k1=row.K1L/row.L,
                                   tilt=row.TILT)

    def make_sextupole(self, row):
        return elements.Sextupole(eid=row.NAME,
                                  l=row.L,
                                  k2=row.K2L/row.L,
                                  tilt=row.TILT)

    def make_octupole(self, row):
        return elements.Octupole(eid=row.NAME,
                                 l=row.L,
                                 k3=row.K3L/row.L,
                                 tilt=row.TILT)

    def make_vkicker(self, row):
        return elements.Vcor(eid=row.NAME,
                             l=row.L,
                             angle=row.VKICK)
    def make_hkicker(self, row):
        return elements.Hcor(eid=row.NAME,
                             l=row.L,
                             angle=row.HKICK)
