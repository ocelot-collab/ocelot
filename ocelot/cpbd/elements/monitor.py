from ocelot.cpbd.elements.element import Element


# to mark locations of bpms and other diagnostics
class Monitor(Element):

    def __init__(self, l=0.0, eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.x_ref = 0.
        self.y_ref = 0.
        self.x = 0.
        self.y = 0.

    def __str__(self):
        s = 'Monitor : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l =%8.4f m\n' % self.l
        return s
