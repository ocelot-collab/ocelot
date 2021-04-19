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
        s = 'Drift : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l =%8.4f m\n' % self.l
        return s
