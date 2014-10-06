___author__ = 'Sergey Tomin'

from numpy import sqrt

# alternative way to get constants
#import scipy.constants as const
#speed_of_light = const.codata.value('speed of light in vacuum')
#m_e_MeV =  const.codata.value('electron mass energy equivalent in MeV')
#pi = const.pi
#h_eV_s = const.codata.value('Planck constant in eV s')
#ro_e = const.codata.value("classical electron radius")

speed_of_light = 299792458.0 #m/s
m_e_MeV = 0.510998928        # MeV

m_e_eV = m_e_MeV*1e+6
m_e_GeV = m_e_MeV*1.e-3      # GeV
pi = 3.14159265359
h_eV_s = 4.135667516e-15    #eV s
ro_e = 2.8179403267e-15
Cgamma = 4.*pi/3.*ro_e/m_e_MeV**3
Cq = 55./(32.*sqrt(3)*2*pi)*h_eV_s*speed_of_light/m_e_eV


Hz = 1.e-9
KHz = 1.e-6
MHz = 1.e-3
GHz = 1.0
THz = 1.e3

V = 1.-9
KV = 1.e-6
MV = 1.e-3
GV = 1.0

eV = 1.e-9
KeV = 1.e-6
MeV = 1.e-3
GeV = 1.0
TeV = 1.e3
