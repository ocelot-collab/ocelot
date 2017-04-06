__author__ = 'Sergey Tomin'
from ocelot.common.globals import *
from desy.flash.lattices.lattice_rf_red import *
import numpy as np

def find_parameters_RF1(E, Es, Ess, Esss):
    '''
    determine RF amplitude and phase from 'knob' values. ACC1 + ACC39
    '''
    f = 1.3e9
    k=2*pi*f/speed_of_light
    n=3 # harmonic cavity
    A0 = E
    B0 = Es
    C0 = Ess
    D0 = Esss
    X1 = (C0+A0*k**2*n**2)/(k**2*(-1+n**2))
    Y1 = -((D0+B0*k**2*n**2)/(k**3*(-1+n**2)))
    X3 = -((C0+A0*k**2)/(k**2*(-1+n**2)))
    Y3 = (D0+B0*k**2)/(k**3*(-n+n**3))
    V1 = np.sqrt(X1**2 + Y1**2)
    f1 = np.arctan(Y1/X1) + pi/2*(1 - np.sign(X1))
    V3 = np.sqrt(X3**2 + Y3**2)
    f3 = np.arctan(Y3/X3) + pi/2*(1 - np.sign(X3))
    phi11=f1*180./pi
    phi13=f3*180./pi-180.
    return (V1, phi11, V3, phi13)

def find_knob_val_RF1(V1, phi1, V39, phi39):
    '''
    determine RF knob values form RF phase/ampl parameters
    '''
    f = 1.3e9
    k=2*pi*f/speed_of_light
    s = 0.0
    E = V1*np.cos(k*s + phi1 * pi / 180.) - V39*np.cos(3*k*s + phi39 * pi / 180.)
    Ep = -k*V1*np.sin(k*s + phi1* pi / 180.) + 3*k*V39*np.sin(3*k*s + phi39 * pi / 180.)
    Epp = -k**2*V1*np.cos(k*s + phi1 * pi / 180.) + 9*k**2*V39*np.cos(3*k*s + phi39 * pi / 180.)
    Eppp = k**3*V1*np.sin(k*s + phi1 * pi / 180.) - 27*k**3* V39*np.sin(3*k*s + phi39 * pi / 180.)

    return (E, Ep, Epp, Eppp)

def find_parameters_RF2(E, Es):
    '''
    determine RF amplitude and phase from 'knob' values. ACC2
    '''
    f = 1.3e9
    k=2*pi*f/speed_of_light
    X1 = E
    Y1 = -Es/k
    V2 = np.sqrt(X1**2 + Y1**2)
    f2 = np.arctan(Y1/X1) + pi/2*(1 - np.sign(X1))
    phi2=f2*180./pi
    return (V2, phi2)

def find_knob_val_RF2(V2, phi2):
    '''
    determine RF knob values form RF phase/ampl parameters
    '''
    f = 1.3e9
    k=2*pi*f/speed_of_light
    s = 0.0
    E = V2*np.cos(k*s + phi2 * pi / 180.)
    Ep = -k*V2*np.sin(k*s + phi2* pi / 180.)
    return (E, Ep)


def rf_cavities(lat):
    big_group = {'name': 'cavity', 'type': 'group', 'children':[]}
    for elem in lat.sequence:
        if elem.__class__ in [Cavity]:
            name = elem.id.split('.')[1]
            if name in [p['name'] for p in big_group['children']]:
                continue
            devs_in_group = []
            tmp = {"name": None, "type":"bool", "value": False}
            tmp["name"] = name + '/SP.AMPL'
            devs_in_group.append(tmp)

            tmp = {"name": None, "type":"bool", "value": False}
            tmp["name"] = name + '/SP.PHASE'
            devs_in_group.append(tmp)

            loc_group = {'name': name, 'type': 'group', 'expanded': False, 'children': devs_in_group}
            big_group['children'].append(loc_group)
    return big_group



def elements(lat, section, groupname='quadrupoles', elemtypes=[Quadrupole]):
    name_section =section["first"] #"GUN-ACC39"
    big_group = {'name': groupname, 'type': 'group', 'children':[]}
    loc_group = {'name': name_section, 'type': 'group', 'expanded': False, 'children':[]}
    devs_in_group = []
    for elem in lat.sequence:
        if elem.__class__ == Marker and elem.id in section.keys():
            name_section = section[elem.id]
            loc_group['children'] = devs_in_group
            big_group['children'].append(loc_group)
            loc_group = {'name': name_section, 'type': 'group', 'expanded': False,'children':[]}
            devs_in_group = []
        elif elem.__class__ in elemtypes:
            if elem.id in [p['name'] for p in devs_in_group]:
                continue
            tmp = {"name": None, "type":"bool", "value": False}
            tmp["name"] = elem.id
            devs_in_group.append(tmp)
    loc_group['children'] = devs_in_group
    big_group['children'].append(loc_group)
    return big_group


def create_group(groupname="group", devnames=['dev1',]):
    group = {'name': groupname, 'type': 'group', 'children':[]}
    for name in devnames:
        group['children'].append({"name": name, "type":"bool", "value": False})
    return group


def generate_tree_params(lat):
    devices = []

    section = {"first": "GUN-ACC39", 'STARTUBC2': "BC2", 'STARTACC2': "ACC23", 'STARTUBC3': "BC3", 'STARTACC4': "ACC4-7", 'ENDACC7':"DOGLEG",
               'STARTSMATCH1': "MATCH", 'STARTUND': "UND"}
    quads = elements(lat, section, groupname='quadrupoles', elemtypes=[Quadrupole])
    cors = elements(lat, section, groupname='correctors', elemtypes=[Hcor, Vcor])
    devices.append(quads)
    devices.append(cors)

    cavs = rf_cavities(lat)
    devices.append(cavs)
    devices.append(create_group(groupname="solenoid", devnames=['1GUN']))
    devices.append(create_group(groupname="dipols", devnames=['D1ECOL']))
    #devices.append(create_group(groupname="SumvoltACC139", devnames=['SUMVOLTAGE.CURVATURE','SUMVOLTAGE.THIRDDERIVATIVE']))
    #devices.append(knobs())
    return devices



def generate_tree_params_old(lat):
    devices = []
    cors = {'name': 'correctors', 'type': 'group', 'children':[]}

    quads = {'name': 'quadrupoles', 'type': 'group', 'children':[]}

    section = {'STARTUBC2':"BC2", 'STARTACC2':"ACC23", 'STARTUBC3':"BC3", 'STARTACC4':"ACC4-7", 'ENDACC7':"DOGLEG",
               'STARTSMATCH1':"MATCH", 'STARTUND':"UND"}

    name_seq = "GUN-ACC39"
    sub_cor_seq = {'name': name_seq, 'type': 'group','expanded': False, 'children':[]}
    sub_cor_chld = []

    sub_quad_seq = {'name': name_seq, 'type': 'group','expanded': False, 'children':[]}
    sub_quad_chld = []

    for elem in lat.sequence:
        if elem.__class__ == Marker and elem.id in section.keys():

            name_seq = section[elem.id]
            #correctors
            sub_cor_seq['children']= sub_cor_chld
            cors['children'].append(sub_cor_seq)
            sub_cor_seq = {'name': name_seq, 'type': 'group','expanded': False,'children':[]}
            sub_cor_chld = []

            #quadrupoles
            sub_quad_seq['children'] = sub_quad_chld
            quads['children'].append(sub_quad_seq)
            sub_quad_seq = {'name': name_seq, 'type': 'group','expanded': False,'children':[]}
            sub_quad_chld = []

        if elem.__class__ in [Hcor, Vcor]:
            tmp = {"name": None, "type":"bool", "value": False}
            tmp["name"] = elem.id
            sub_cor_chld.append(tmp)

        elif elem.__class__ in [Quadrupole]:
            if elem.id in [p['name'] for p in sub_quad_chld]:
                continue
            tmp = {"name": None, "type":"bool", "value": False}
            tmp["name"] = elem.id
            sub_quad_chld.append(tmp)

    sub_cor_seq['children'] = sub_cor_chld
    cors['children'].append(sub_cor_seq)

    sub_quad_seq['children'] = sub_quad_chld
    quads['children'].append(sub_quad_seq)
    devices.append(cors)
    devices.append(quads)
    return devices
