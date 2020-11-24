from ocelot.cpbd.transformations.optics import t_nnn, SecondOrderMult
from ocelot.cpbd.transformations.corrector import CorrectorTM


class VCorrectorTM(CorrectorTM):
    @classmethod
    def create_from_element(cls, element, params=None):
        t_mat_z_e = lambda z, energy: t_nnn(z, 0, 0, 0, energy)
        tm = cls(angle_x=0, angle_y=element.angle, r_z_no_tilt=element.create_r_matrix(), t_mat_z_e=t_mat_z_e)
        tm.multiplication = SecondOrderMult().tmat_multip
        return tm