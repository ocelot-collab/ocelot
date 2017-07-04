___author__ = 'Sergey Tomin'

from numpy import sqrt

# alternative way to get constants
#import scipy.constants as const
#speed_of_light = const.codata.value('speed of light in vacuum')
#m_e_MeV =  const.codata.value('electron mass energy equivalent in MeV')
#pi = const.pi
#h_eV_s = const.codata.value('Planck constant in eV s')
#ro_e = const.codata.value("classical electron radius")

pi = 3.141592653589793
speed_of_light = 299792458.0 #m/s
q_e = 1.6021766208e-19       # C - Elementary charge
m_e_kg = 9.10938215e-31      # kg
m_e_eV = m_e_kg * speed_of_light**2 / q_e # eV (510998.8671)
m_e_MeV = m_e_eV / 1e+6        # MeV (0.510998928)
m_e_GeV = m_e_eV / 1e+9      # GeV


mu_0 = 4*pi*1e-7
epsilon_0 = 8.854187817620e-12 # F/m

h_eV_s = 4.135667516e-15    #eV s
hr_eV_s = h_eV_s/2./pi
ro_e = 2.8179403267e-15     # m

Cgamma = 4.*pi/3.*ro_e/m_e_MeV**3
Cq = 55./(32.*sqrt(3)*2*pi)*h_eV_s*speed_of_light/m_e_eV

Z0 = 1./(speed_of_light*epsilon_0)  # Ohm - impedance of free space

"""
def lambda2eV(Lambda):
    Eph = h_eV_s*speed_of_light/Lambda
    return Eph

def eV2lambda(Ephoton):
    Lambda = h_eV_s*speed_of_light/Ephoton
    return Lambda

def Ephoton2K(Eph, lu = 0.04, Eeb = 14):
    gamma = Eeb/m_e_GeV
    K = sqrt(4.*gamma*gamma/lu*h_eV_s*speed_of_light/Eph - 2)
    return K

def K2Ephoton(K, lu = 0.04, E=14):
    gamma = E/m_e_GeV
    Eph = 4.*gamma*gamma*h_eV_s*speed_of_light/((K*K + 2)*lu)
    return Eph

def K2Lambda(K, lu = 0.04, E=14):
    gamma = E/m_e_GeV
    Eph = 4.*gamma*gamma*h_eV_s*speed_of_light/((K*K + 2)*lu)
    Lambda = eV2lambda(Eph)
    return Lambda

def field2K(field, lu = 0.04):
    K = field*lu*speed_of_light/(m_e_eV*2.*pi)
    return K
print (field2K(0.65, lu = 0.007) )
#print  0.66*0.007*1.6e-19/(9.1e-31*speed_of_light*2.*pi)
def K2field(K, lu = 0.04):
    field = K*m_e_eV*2.*pi/(lu*speed_of_light)
    return field

def field2Ephoton(field, lu = 0.04, E=14):
    gamma = E/m_e_GeV
    K = field2K(field, lu)
    l = lu/(2.*gamma*gamma)*(1.+K*K/2.)
    Eph = h_eV_s*speed_of_light/l
    return Eph

def Ephoton2field(energy, lu = 0.04, Eeb = 14):
    K = Ephoton2K(energy, lu, Eeb)
    field = K*2.*pi*m_e_eV/(lu*speed_of_light)
    return field

"""



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
