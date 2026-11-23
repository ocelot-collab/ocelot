from ocelot.cpbd.transformations.transfer_map import TransferMap, TMTypes
from ocelot.cpbd.elements.element import Element


class RungeKuttaGlobalTM(TransferMap):
    """Track particles between fixed Cartesian Z planes in a magnetic field.

    Use this method to obtain positions in the field frame, for example when
    comparing an orbit to a measured field map. A reference electron bends
    in a dipole even if all its entrance coordinates are zero. ParticleArray
    stores mechanical transverse momenta divided by the reference momentum.
    Arrival time is tau=ct-Z/beta0 in metres, including the entrance delay.

    Use RungeKuttaOcelotTM to follow the beam in the moving reference frame
    through a bend. Laboratory placement of successive elements requires
    their survey geometry.
    """

    def __init__(self, create_tm_param_func, delta_e_func, tm_type: TMTypes, length: float, delta_length: float = 0.0, **params) -> None:
        super().__init__(create_tm_param_func, delta_e_func, tm_type, length, delta_length)
        s_start = params.get('s_start')
        self.s_start = s_start if s_start else 0
        npoints = params.get('npoints')
        self.npoints = npoints if npoints is not None else 200
        self.long_dynamics = True

    @classmethod
    def from_element(cls, element: Element, tm_type: TMTypes = TMTypes.MAIN, delta_l=None, **params):
        """Bind the field callback for the requested element body or edge."""
        return cls.create(entrance_tm_params_func=element.create_runge_kutta_entrance_params if element.has_edge else None,
                          delta_e_func=element.create_delta_e,
                          main_tm_params_func=element.create_runge_kutta_main_params,
                          exit_tm_params_func=element.create_runge_kutta_exit_params if element.has_edge else None,
                          tm_type=tm_type, length=element.l, delta_length=delta_l, **params)

    def get_params(self, energy):
        return self.create_tm_param_func(energy)

    def map_function(self, X, energy: float):
        from ocelot.cpbd.high_order import rk_field

        params = self.get_params(energy)
        length = self.delta_length if self.delta_length is not None else self.length
        return rk_field(X, self.s_start, self.s_start + length, self.npoints, energy, params.mag_field, self.long_dynamics)


class RungeKuttaTM(RungeKuttaGlobalTM):
    """
    Legacy name for :class:`RungeKuttaGlobalTM`.

    Kept for backwards compatibility with existing Ocelot examples and user
    lattices.  New beam tracking code that expects Ocelot-relative coordinates
    should use ``RungeKuttaOcelotTM``.
    """
    pass


class RungeKuttaOcelotTM(RungeKuttaGlobalTM):
    """
    Runge-Kutta transformation that returns Ocelot-relative beam coordinates.

    The field integration is Cartesian, while input/output coordinates belong
    to normal planes of one reference defined at the element entrance. Slices
    retain that reference and its Bishop axes. Transverse momenta are mechanical
    p/p0; no vector-potential terms are added. Arrival time is tau=ct-s/beta0.
    npoints bounds the RK step by full_element_length/(npoints-1).

    reference="field" follows the native field reference (the default).
    reference="nominal" uses the surveyed element frames to measure orbit errors.
    """
    reference = "field"

    def __init__(self, *args, reference="field", nominal_face=None, **params):
        super().__init__(*args, **params)
        if reference not in ("field", "nominal"):
            raise ValueError("RK reference must be 'field' or 'nominal'")
        if reference == "nominal" and nominal_face is None:
            raise ValueError("Nominal RK reference requires the element's survey geometry")
        self.reference = reference
        self._nominal_face = nominal_face

    @classmethod
    def from_element(cls, element, tm_type=TMTypes.MAIN, delta_l=None, **params):
        if params.get("reference") == "nominal":
            params["nominal_face"] = element._rk_nominal_face
        return super().from_element(element, tm_type, delta_l, **params)

    def map_function(self, X, energy: float):
        from ocelot.cpbd.high_order import rk_field_to_ocelot

        params = self.get_params(energy)
        length = self.delta_length if self.delta_length is not None else self.length
        alignment = getattr(params, "alignment", None)
        field = params.design_field if alignment else params.mag_field
        use_local_offsets = self.reference == "nominal" and getattr(params, "centered_field", None) is not None
        if use_local_offsets:
            field = params.centered_field
            alignment = params.nominal_alignment
        return rk_field_to_ocelot(X, self.s_start, self.s_start + length, self.npoints,
                                  energy, field,
                                  self.long_dynamics, reference_length=self.length, alignment=alignment,
                                  local_offsets=use_local_offsets,
                                  reference_face=self._nominal_face if self.reference == "nominal" else None)
