from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.unknown_element_atom import UnknownElementAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap


class UnknownElement(OpticElement):
    """
    l - length of lens in [m]
    """
    def __init__(self, l=0, kick=0, xsize=0, ysize=0, volt=0, lag=0, harmon=0, refer=0, vkick=0, hkick=0, eid=None,
                 tm=TransferMap, **params):
        super().__init__(UnknownElementAtom(l=l, kick=kick, xsize=xsize, ysize=ysize, volt=volt, lag=lag, harmon=harmon,
                                            refer=refer, vkick=vkick, hkick=hkick, eid=eid),
                         tm=tm, default_tm=TransferMap, **params)