from ocelot.cpbd.tm_params.tm_params import TMParams


class RungeKuttaParams(TMParams):
    """
    Field callback container for ``RungeKuttaTM``.

    ``mag_field`` is a callable returning ``(Bx, By, Bz)`` for the current
    coordinates, so the transformation can integrate through the field.
    """

    def __init__(self, mag_field, design_field=None, alignment=None) -> None:
        super().__init__()
        self.mag_field = mag_field
        self.design_field = mag_field if design_field is None else design_field
        self.alignment = alignment
    
