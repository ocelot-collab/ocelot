from ocelot.cpbd.elements.magnet import Magnet


class DriftAtom(Magnet):
    """
    drift - free space
    l - length of drift in [m]
    """

    def __init__(self, l=0., eid=None):
        super().__init__(eid)
        self.l = l

    def __str__(self):
        s = 'Drift('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s