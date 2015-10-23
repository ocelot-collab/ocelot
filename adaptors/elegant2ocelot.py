'''
Created on 05.04.2013
@author: zagor
'''
import csv
from ocelot.cpbd.elements import *
from math import *

def read_lattice_elegant(file_flo, file_par):
    f=open(file_flo,'rb')
    data_flo=csv.reader(f,delimiter='\t')
    data_flo=[row for row in data_flo]
    f.close()
    f=open(file_par,'rb')
    data_par=csv.reader(f, delimiter='\t')
    data_par=[row for row in data_par]
    f.close()
    lattice=[]
    n_flo=len(data_flo)
    i_s=0;# i_X=1; i_Y=2; 
    i_Z=3;#i_theta=4;i_phi=5;    i_psi=6;    
    i_ElementName=7;# i_ElementOccurence=8; 
    i_ElementType=9; 
    for i in range(8,n_flo):
        v=data_flo[i] 
        sname=v[i_ElementName]
        stype=v[i_ElementType];
        sname=sname.replace('-C','')
        if stype=='QUAD':
            quad = Quadrupole(id=sname)
            quad.s=eval(v[i_s])
            quad.z=eval(v[i_Z])
            lattice=lattice+[quad]
        elif stype=='SEXT':
            sext = Sextupole(id=sname)
            sext.s=eval(v[i_s])
            sext.z=eval(v[i_Z])
            lattice=lattice+[sext]
        elif (stype=='CSBEND') or (stype == 'CSRCSBEND') or (stype == 'SBEN'):
            sben = Bend(l=1,id=sname)
            sben.s=eval(v[i_s])
            sben.z=eval(v[i_Z])
            lattice=lattice+[sben]
        elif (stype=='CRBEND') or (stype=='CSRCRBEND'):
            rben= RBend (id=sname)
            rben.s=eval(v[i_s])
            rben.z=eval(v[i_Z])
            lattice=lattice+[rben]
        elif stype=='WIGGLER':
            undulator = Undulator(id=sname, lperiod=0, nperiods=0, Kx=0)
            undulator.s=eval(v[i_s])
            undulator.z=eval(v[i_Z])
            lattice=lattice+[undulator]   
        elif stype=='RFCA':
            cavity = Cavity(id=sname, l=0)
            cavity.s=eval(v[i_s])
            cavity.z=eval(v[i_Z])
            lattice=lattice+[cavity]   
        elif stype=='HKICK':
            hcor = Hcor(id=sname)
            hcor.s=eval(v[i_s])
            hcor.z=eval(v[i_Z])
            lattice=lattice+[hcor]
        elif stype=='VKICK':
            vcor = Vcor(id=sname)
            vcor.s=eval(v[i_s])
            vcor.z=eval(v[i_Z])
            lattice=lattice+[vcor]
        elif stype=='MONI':
            monitor = Monitor(id=sname)
            monitor.s=eval(v[i_s])
            monitor.z=eval(v[i_Z])
            lattice=lattice+[monitor]
            
    n_par=len(data_par)
    pos=6
    for elem in lattice:
        name0=elem.id
        name1=data_par[pos][0]
        while (pos<n_par-1) and (name0!=name1):
            pos=pos+1
            name1=data_par[pos][0]
            stype=data_par[pos][3]
        if pos>n_par-2:
            break
        elem.l=eval(data_par[pos][2])
        if elem.type=="quadrupole":
            elem.k1=eval(data_par[pos+1][2])
        elif elem.type=="sextupole":
            elem.k2=eval(data_par[pos+1][2])
        elif elem.type=="bend":
            elem.angle=eval(data_par[pos+1][2])
            if (stype=='CRBEND') or (stype=='CSRCRBEND') or (stype=='CSRCSBEND'):
                elem.e1=eval(data_par[pos+10][2])
                elem.e2=eval(data_par[pos+11][2])
                elem.tilt=eval(data_par[pos+12][2])
            if stype=='SBEN':
                elem.e1=eval(data_par[pos+3][2])
                elem.e2=eval(data_par[pos+4][2])
                elem.tilt=eval(data_par[pos+5][2])
        elif elem.type=="rbend":
            elem.angle=eval(data_par[pos+1][2])
            elem.tilt=eval(data_par[pos+12][2])
        elif elem.type=="undulator":
            elem.l=eval(data_par[pos][2])
            elem.nperiods=eval(data_par[pos+7][2])/2
            elem.lperiod=elem.l/elem.nperiods
            elem.Kx=eval(data_par[pos+2][2])
            elem.Ky=0
        elif elem.type=="cavity":
            elem.l=eval(data_par[pos][2])
            elem.v=eval(data_par[pos+1][2])
            elem.phi=eval(data_par[pos+2][2])-90
            elem.f=eval(data_par[pos+3][2])
            elem.delta_e=elem.v*1e-9*cos(elem.phi*pi/180)        
        elif elem.type=="hcor":
            elem.l=eval(data_par[pos][2])
        elif elem.type=="vcor":
            elem.l=eval(data_par[pos][2])
        elif elem.type=="monitor":
            elem.l=eval(data_par[pos][2]) 
        elem.id=elem.id.replace('.','_')              
    return lattice

def insert_drifts(z_start, z_stop, lat_def):
    lattice= []
    s0=z_start
    for elem in lat_def:
        if elem.z>z_start and elem.z<z_stop:
            if s0==z_start:
                ds=elem.z-elem.l-s0
            else:
                ds=elem.s-elem.l-s0
            if ds>0:
                str_id=str(ds).replace('.','_')
                lattice=lattice+[Drift(l=ds,id='DRIFT_'+str_id)]
            lattice=lattice+[elem]
            print elem.id,elem.type,elem.z, elem.s, ds
            s0=elem.s
        elif elem.z>z_stop:
            break
    ds=z_stop-elem.z
    if ds>0:
        str_id=str(ds).replace('.','_')
        lattice=lattice+[Drift(l=ds,id='DRIFT'+str_id)]
    return lattice

    
     


    

