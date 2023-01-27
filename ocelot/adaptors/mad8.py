__author__ = "Stuart Walker"

import contextlib
import copy

import numpy as np
import pandas as pd

from ocelot.cpbd.beam import Twiss
import ocelot.cpbd.elements as elements
from ocelot.common.globals import m_e_GeV, speed_of_light

try:
    import pand8
except ImportError:
    pass


TWISS_MAD8_TO_OCELOT_ELEMENT_NAMES = {
    "DRIF": "Drift",
    "HKIC": "Hcor",
    # 'KICK': "Kicker", Must be split into two kickers
    "MARK": "Marker",
    "MONI": "Monitor",
    "QUAD": "Quadrupole",
    "RBEN": "RBend",
    "SBEN": "SBend",
    "VKIC": "Vcor",
}

TWISS_MAD8_TO_OCELOT_KEYS = {
    "NAME": "eid",
    "H1": "h_pole1",
    "H2": "h_pole2",
}


def twiss_to_sequence_with_optics(twiss_file):
    try:
        twiss_df = pand8.read(twiss_file)
    except TypeError:
        # then assume it's already a pandas DataFrame
        twiss_df = twiss_file

    twiss_df = convert_twiss_kickers_to_hvkicker_pairs(twiss_df)
    twiss0 = twiss_to_ocelot_twiss0(twiss_df)
    sequence = twiss_to_ocelot_sequence(twiss_df)
    return sequence, twiss0


def convert_twiss_kickers_to_hvkicker_pairs(twiss):
    """Replace KICK elements in the MAD8 Twiss DataFrame with a pair of equivalent HKIC
    and VKIC elements.

    """

    kick_indices = twiss.index[twiss["KEYWORD"] == "KICK"]
    while len(kick_indices):
        this_index = kick_indices[0]

        kicker = twiss.iloc[this_index]
        before = twiss.iloc[:this_index]
        after = twiss.iloc[this_index + 1 :]

        # The two kickers to replace the original one.
        hkicker = kicker.copy(deep=True)
        vkicker = kicker.copy(deep=True)
        hkicker["NAME"] += "_H"
        vkicker["NAME"] += "_V"
        hkicker["KEYWORD"] = "HKIC"
        vkicker["KEYWORD"] = "VKIC"
        hkicker["VKICK"] = 0.0
        hkicker["ANGLE"] = hkicker["HKICK"]
        vkicker["HKICK"] = 0.0
        vkicker["ANGLE"] = vkicker["VKICK"]
        hkicker["L"] /= 2.0
        vkicker["L"] /= 2.0

        both_new_kickers = pd.DataFrame([hkicker, vkicker])

        twiss = pd.concat((before, both_new_kickers, after)).reset_index(drop=True)
        kick_indices = twiss.index[twiss["KEYWORD"] == "KICK"]

    return twiss


def twiss_to_ocelot_twiss0(twiss_df, exn=1e-6, eyn=1e-6):
    out = Twiss()

    # For some reason E is 0 for the first row, so use the next row to get E...
    out.E = twiss_df.iloc[1].E
    rel_gamma = out.E / m_e_GeV
    rel_beta = np.sqrt(1 - rel_gamma**-2)

    out.emit_xn = exn  # out.emit_x * rel_gamma
    out.emit_yn = eyn  # out.emit_y * rel_gamma
    out.emit_x = out.emit_x / rel_gamma / rel_beta
    out.emit_y = out.emit_y / rel_gamma / rel_beta

    first = twiss_df.iloc[0]
    out.alpha_x = first["ALFX"]
    out.alpha_y = first["ALFY"]
    out.beta_x = first["BETX"]
    out.beta_y = first["BETY"]
    out.Dx = first["DX"]
    out.Dy = first["DY"]
    out.Dxp = first["DPX"]
    out.Dyp = first["DPY"]
    out.x = first["X"]
    out.y = first["Y"]
    out.xp = first["PX"]
    out.yp = first["PY"]

    return out


def twiss_to_ocelot_sequence(twiss_df):
    converter = TwissToOCELOTConverter()
    for row in twiss_df.itertuples():
        if not row.KEYWORD and row.NAME == "INITIAL":
            continue
        yield converter.dispatch(row)


class UnsupportedElementError(Exception):
    pass


class TwissToOCELOTConverter:
    MAD8_TWISS_NAMES_TO_OCELOT = {
        "DRIF": "Drift",
        "HKIC": "Hcor",
        "VKIC": "Vcor",
        "MARK": "Marker",
        "MONI": "Monitor",
        "QUAD": "Quadrupole",
        "RBEN": "RBend",
        "SBEN": "SBend",
    }

    def dispatch(self, row):
        etype = row.KEYWORD
        name = row.NAME
        try:
            ocelot_class_name = self.MAD8_TWISS_NAMES_TO_OCELOT[etype]
        except:
            raise UnsupportedElementError(
                f"Unable to map {etype} to an OCELOT class name"
            )

        try:
            return getattr(self, f"make_{ocelot_class_name}")(name, row)
        except AttributeError:
            if etype == "KICK":
                raise UnsupportedElementError(
                    "Kicker elements must be split into HKICKER and VKICKER pairs."
                )

            raise UnsupportedElementError(
                f"Missing construction method for OCELOT"
                f" element type: {ocelot_class_name}."
            )

    def make_Drift(self, name, attributes):
        return elements.Drift(eid=name, l=attributes.L)

    def make_Hcor(self, name, attributes):
        return elements.Hcor(eid=name, angle=attributes.HKICK, l=attributes.L)

    def make_Vcor(self, name, attributes):
        return elements.Vcor(eid=name, angle=attributes.ANGLE, l=attributes.L)

    def make_Marker(self, name, attributes):
        return elements.Marker(eid=name)

    def make_Monitor(self, name, attributes):
        return elements.Monitor(eid=name, l=attributes.L)

    def make_Quadrupole(self, name, attributes):
        at = attributes
        return elements.Quadrupole(eid=name, l=at.L, k1=at.K1, tilt=at.TILT)

    def make_RBend(self, name, attributes):
        at = attributes
        return elements.RBend(
            eid=name,
            l=at.L,
            angle=at.ANGLE,
            k1=at.K1,
            e1=at.E1,
            e2=at.E2,
            tilt=at.TILT,
            k2=at.K2,
            h_pole1=at.H1,
            h_pole2=at.H2,
        )

    def make_SBend(self, name, attributes):
        at = attributes
        return elements.SBend(
            eid=name,
            l=at.L,
            angle=at.ANGLE,
            k1=at.K1,
            e1=at.E1,
            e2=at.E2,
            tilt=at.TILT,
            k2=at.K2,
            h_pole1=at.H1,
            h_pole2=at.H2,
        )
        return elements.SBend(eid=name)
