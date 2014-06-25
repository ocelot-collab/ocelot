'''
I/O operations
'''
from numpy import shape
import h5py
from ocelot.cpbd.beam import Beam
from ocelot.common.screen import Screen

class XIO:
    def __init__(self, filename, mode='w'):
        self.filename = filename
        self.file = h5py.File(filename,mode)

    def close(self):
        self.file.close()


    def write_beam(self, beam):

        gi = self.file.create_group('Beam')
        gi.attrs['E'] = beam.E
        gi.attrs['I'] = beam.I
        gi.attrs['sigma_E'] = beam.sigma_E
        gi.attrs['emit_x'] = beam.emit_x
        gi.attrs['emit_y'] = beam.emit_y
        gi.attrs['beta_x'] = beam.beta_x
        gi.attrs['beta_y'] = beam.beta_y
        gi.attrs['alpha_x'] = beam.alpha_x
        gi.attrs['alpha_y'] = beam.alpha_y

    def read_beam(self):
        if not( 'Beam' in self.file.keys()):
            return None
        gi = self.file['Beam']
        beam = Beam()

        beam.E = gi.attrs['E']
        beam.I = gi.attrs['I']
        beam.sigma_E = gi.attrs['sigma_E']
        beam.emit_x = gi.attrs['emit_x']
        beam.emit_y = gi.attrs['emit_y']
        beam.beta_x = gi.attrs['beta_x']
        beam.beta_y = gi.attrs['beta_y']
        beam.alpha_x = gi.attrs['alpha_x']
        beam.alpha_y = gi.attrs['alpha_y']

        return beam

    def write_screen(self, screen):

        gi = self.file.create_group('Screen')
        gi.attrs['z'] = screen.z
        gi.attrs['size_x'] = screen.size_x
        gi.attrs['size_y'] = screen.size_y
        gi.attrs['nx'] = screen.nx
        gi.attrs['ny'] = screen.ny
        gi.attrs['start_energy'] = screen.start_energy
        gi.attrs['end_energy'] = screen.end_energy
        gi.attrs['num_energy'] = screen.num_energy

    def read_screen(self):
        if not('Screen' in self.file.keys()):
            return None
        gi = self.file['Screen']
        screen = Screen()
        screen.z = gi.attrs['z']
        screen.size_x = gi.attrs['size_x']
        screen.size_y = gi.attrs['size_y']
        screen.nx = gi.attrs['nx']
        screen.ny = gi.attrs['ny']
        screen.start_energy = gi.attrs['start_energy']
        screen.end_energy = gi.attrs['end_energy']
        screen.num_energy = gi.attrs['num_energy']

        return screen

    def write_intensity(self, i, idx = 0):

        gi = self.file.create_group('Intensities')
        ds = gi.create_dataset('I_' + str(idx), (shape(i)), 'f', chunks=True)
        ds[:] = i
        ds.attrs['energy_unit'] =  'ev'

    def read_intensity(self, idx = 0):
        if not('Intensities' in self.file.keys()):
            return None
        gi = self.file['Intensities']
        intens = gi['I_0'].value
        return intens

    def write_trajectory(self, tr, idx = 0):

        gt = self.file.create_group('Trajectory')
        gt2 = gt.create_group('T_'+ str(idx))
        ds_x = gt2.create_dataset('T_X', (shape(tr.arX)), 'f', chunks=True)
        ds_x[:] = tr.arX
        ds_xp = gt2.create_dataset('T_Xp', (shape(tr.arXp)), 'f', chunks=True)
        ds_xp[:] = tr.arXp
        ds_y = gt2.create_dataset('T_Y', (shape(tr.arY)), 'f', chunks=True)
        ds_y[:] = tr.arY
        ds_yp = gt2.create_dataset('T_Yp', (shape(tr.arYp)), 'f', chunks=True)
        ds_yp[:] = tr.arYp
        ds_z = gt2.create_dataset('T_Z', (shape(tr.arZ)), 'f', chunks=True)
        ds_z[:] = tr.arZ
        ds_zp = gt2.create_dataset('T_Zp', (shape(tr.arZp)), 'f', chunks=True)
        ds_zp[:] = tr.arZp

        gt.attrs['npoints'] = len(tr.arX)

    def read_trajectory(self, idx = 0):
        if not('Trajectory' in self.file.keys()):
            return None
        gi = self.file['Trajectory']
        gt = gi['T_0']
        x = gt['T_X'].value
        xp = gt['T_Xp'].value
        y = gt['T_Y'].value
        yp = gt['T_Yp'].value
        z = gt['T_Z'].value
        zp = gt['T_Zp'].value

        return (x,xp,y,yp,z,zp)


class Dump:
    def __init__(self):
        self.readme = ''
        self.index = {}

    def dump(self, xio):
        if self.index.has_key('beam'):
            xio.write_beam(self.index['beam'])

        if self.index.has_key('screen'):
            xio.write_screen(self.index['screen'])

        if self.index.has_key('intensity'):
            xio.write_intensity(self.index['intensity'])

        if self.index.has_key('trajectory'):
            xio.write_trajectory(self.index['trajectory'])


        xio.close()

    def read(self, xio):
        self.index['beam'] = xio.read_beam()
        self.index['screen'] = xio.read_screen()
        self.index['intensity'] = xio.read_intensity()
        self.index['trajectory'] = xio.read_trajectory()
