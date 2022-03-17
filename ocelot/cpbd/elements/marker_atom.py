from ocelot.cpbd.elements.element import Element


class MarkerAtom(Element):
    def __init__(self, eid=None):
        Element.__init__(self, eid)
        self.l = 0.

    def __str__(self):
        s = 'Marker('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s