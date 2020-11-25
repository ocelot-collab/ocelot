from ocelot.cpbd.elements.element import Element


class Marker(Element):
    def __init__(self, eid=None):
        Element.__init__(self, eid)
        self.l = 0.

    def __str__(self):
        s = 'Marker : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l =%8.4f m\n' % self.l
        return s