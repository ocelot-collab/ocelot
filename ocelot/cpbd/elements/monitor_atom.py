from ocelot.cpbd.elements.element import Element


# to mark locations of bpms and other diagnostics
class MonitorAtom(Element):

    def __init__(self, l=0.0, eid=None):
        super().__init__(eid)
        self.l = l
        self.x_ref = 0.
        self.y_ref = 0.
        self.x = 0.
        self.y = 0.

    def __str__(self):
        s = 'Monitor('
        s += 'l=%7.5f, ' % self.l if self.l != 0. else ""
        s += 'eid="' + str(self.id) + '")' if self.id is not None else ")"
        return s