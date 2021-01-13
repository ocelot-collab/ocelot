from ocelot.cpbd.elements.element import Element


class UnknownElement(Element):
    """
    l - length of lens in [m]
    """

    def __init__(self, l=0, kick=0, xsize=0, ysize=0, volt=0, lag=0, harmon=0, refer=0, vkick=0, hkick=0, eid=None):
        Element.__init__(self, eid)
        self.l = l