
from numpy import sqrt

speed_of_light = 299792458.0 #m/s
m_e_MeV = 0.510998928        # MeV
m_e_eV = m_e_MeV*1e+6
m_e_GeV = m_e_MeV*1.e-3      # GeV
h_eV_s = 4.135667516e-15    #eV s
pi = 3.14159265359

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
#print (field2K(0.65, lu = 0.007) )
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


