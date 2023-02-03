from ocelot.cpbd.elements.magnet import Magnet


class UnknownElementAtom(Magnet):
    """
    Unknown element
    """

    def __init__(self, l=0, kick=0, xsize=0, ysize=0, volt=0, lag=0, harmon=0, refer=0, vkick=0, hkick=0, eid=None):
        super().__init__(eid)

    def __str__(self):
        s = 'UnknownElement('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s