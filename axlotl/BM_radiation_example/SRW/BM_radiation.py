# -*- coding: utf-8 -*-
#############################################################################
# SRWLIB Example#13: Simulating emission of Bending Magnet SR
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from uti_plot import * #required for plotting
import time
import numpy as np
import os

print('SRWLIB Python Example # 13:')
print('Simulating emission of Bending Magnet Synchrotron Radiaiton')

#***********Data Folder and File Names
FolderName = 'BM_radiation' #example data sub-folder name
FieldX_fromElectron = 'SRW_Ex_field_from_electron.txt'
FieldY_fromElectron = 'SRW_Ey_field_from_electron.txt'
FieldX_fromBeam = 'SRW_Ex_field_from_beam.txt' #file name for output propagated single-electron SR intensity vs X and Y data
FieldY_fromBeam = 'SRW_Ey_field_from_beam.txt'  #file name for output propagated multi-electron SR intensity vs X and Y data

#***********Bending Magnet
B = 0.10439 #Dipole magnetic field [T]
LeffBM = 5.378 #Magnet length [m]
BM = SRWLMagFldM(B, 1, 'n', LeffBM)
magFldCnt = SRWLMagFldC([BM], [0], [0], [0]) #Container of magnetic field elements and their positions in 3D

#***********Electron Beam
eBeam = SRWLPartBeam()
Iavg = 0.1 #Average current [A]
Energy = 6.08 #Beam Energy [GeV]
En_spread = 0.001 #Energy spread
Hor_emit = 1e-9 # Horizontal emittance [m]
Hor_beta = 28.776 # Horizontal Beta[m]
Hor_alpha = 0 # Horizontal Alpha [rad]
Hor_disp = 0.005796 # Horizontal dispersion function [m]
Hor_der = 0.0013 # Horizontal derivative [rad]
Vrt_emit = 0.01e-9 # Vertical emittance [m]
Vrt_beta = 9.508 # Vertical Beta [m]
Vrt_alpha = 0 # Vertical Alpha [rad]
Vrt_disp = 0 # Vertical dispersion function [m]
Vrt_der = 0 # Vertical derivative [rad]
#1st order statistical moments:
eBeam.partStatMom1.x = 0. #Initial horizontal position of central trajectory [m]
eBeam.partStatMom1.y = 0. #Initial vertical position of central trajectory [m]
eBeam.partStatMom1.z = 0. #Initial longitudinal position of central trajectory [m]
eBeam.partStatMom1.xp = 0. #Initial horizontal angle of central trajectory [rad]
eBeam.partStatMom1.yp = 0. #Initial vertical angle of central trajectory [rad]
eBeam.from_Twiss(Iavg, Energy, En_spread, Hor_emit, Hor_beta, Hor_alpha, Hor_disp, Hor_der, Vrt_emit, Vrt_beta, Vrt_alpha, Vrt_disp, Vrt_der) #Initial Beam parameters via Twiss parameters

#***********Radiation Sampling for the On-Axis SR Spectrum
'''
wfrSp = SRWLWfr() #Wavefront structure (placeholder for data to be calculated)
wfrSp.allocate(500, 1, 1) #Numbers of points vs photon energy, horizontal and vertical positions (the last two will be modified in the process of calculation)
wfrSp.mesh.zStart = 5. #Longitudinal position for initial wavefront [m]

wfrSp.mesh.eStart = 0.1 #Initial photon energy [eV]
wfrSp.mesh.eFin = 10000. #Final photon energy [eV]

wfrSp.mesh.xStart = 0. #Initial horizontal position [m]
wfrSp.mesh.xFin = wfrSp.mesh.xStart #Final horizontal position [m]
wfrSp.mesh.yStart = 0. #Initial vertical position [m]
wfrSp.mesh.yFin = 0. #Final vertical position [m]

wfrSp.partBeam = eBeam #e-beam data is contained inside the wavefront struct
'''

#***********Radiation Sampling for the Initial Wavefront (before first optical element)
wfr = SRWLWfr() #Wavefront structure (placeholder for data to be calculated)
wfr.allocate(1, 101, 101) #Numbers of points vs photon energy, horizontal and vertical positions (the last two will be modified in the process of calculation)

distSrc = 14.2 #Distance from geometrical source point to lens [m]
wfr.mesh.zStart = distSrc #Longitudinal position for initial wavefront [m]

Lambda = 450e-9 # Observed wavelength [m]
c = 2.99e+8 # Speed of light [m/s]
h = 4.135667e-15 # Plank`s constant
Ph_en = (h*c)/(Lambda) # Oserved photons energy [eV]
print('Lambda = ', Lambda*1e+9,'nm')
print('Photon Energy =', Ph_en, 'eV')

wfr.mesh.eStart = Ph_en #Initial photon energy [eV]
wfr.mesh.eFin = wfr.mesh.eStart #Final photon energy [eV]

wfr.mesh.xStart = -0.01 #Initial horizontal position [m]
wfr.mesh.xFin = 0.01 #Final horizontal position [m]
wfr.mesh.yStart = -0.01 #Initial vertical position [m]
wfr.mesh.yFin = 0.01 #Final vertical position [m]

wfr.partBeam = eBeam #e-beam data is contained inside the wavefront struct

#***********BM SR Calculation
#Precision parameters
meth = 2 #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
relPrec = 0.005 #Relative precision
zStartInteg = -5.378/2 #Longitudinal position to start integration (effective if < zEndInteg)
zEndInteg = 5.378/2 #Longitudinal position to finish integration (effective if > zStartInteg)
npTraj = 20000 #Number of points for trajectory calculation 
useTermin = 1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)

print('   Performing initial electric field wavefront calculation ... ', end='')
t0 = time.time()
sampFactNxNyForProp = 0 #Sampling factor for adjusting nx, ny (effective if > 0)
arPrecSR = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, sampFactNxNyForProp]
srwl.CalcElecFieldSR(wfr, 0, magFldCnt, arPrecSR) #Calculating electric field
print('done in', round(time.time() - t0), 's')

print('   Extracting intensity ... ', end='')
t0 = time.time()
mesh0 = deepcopy(wfr.mesh)
arI0 = array('f', [0]*mesh0.nx*mesh0.ny) #"Flat" array to take 2D intensity data (vs X & Y)
arElecEx = array('f', [0]*mesh0.nx*mesh0.ny)
arElecEy = array('f', [0]*mesh0.nx*mesh0.ny)
arBeamEx = array('f', [0]*mesh0.nx*mesh0.ny)
arBeamEy = array('f', [0]*mesh0.nx*mesh0.ny)
srwl.CalcIntFromElecField(arI0, wfr, 6, 0, 3, mesh0.eStart, 0, 0) #Extracting intensity vs horizontal and vertical positions
srwl.CalcIntFromElecField(arElecEx, wfr, 0, 0, 3, mesh0.eStart, 0, 0)
arElecEx1 = np.asarray(arElecEx).reshape((mesh0.nx,mesh0.ny))
srwl.CalcIntFromElecField(arElecEy, wfr, 1, 0, 3, mesh0.eStart, 0, 0)
arElecEy1 = np.asarray(arElecEy).reshape((mesh0.nx,mesh0.ny))
srwl.CalcIntFromElecField(arBeamEx, wfr, 0, 1, 3, mesh0.eStart, 0, 0)
arBeamEx1 = np.asarray(arBeamEx).reshape((mesh0.nx,mesh0.ny))
srwl.CalcIntFromElecField(arBeamEy, wfr, 1, 1, 3, mesh0.eStart, 0, 0)
arBeamEy1 = np.asarray(arBeamEy).reshape((mesh0.nx,mesh0.ny))
print('done in', round(time.time() - t0), 's')

print('   Saving intensity ... ', end='')
t0 = time.time()
np.savetxt(os.path.join(os.getcwd(),FolderName,FieldX_fromElectron), arElecEx1, newline="\n", delimiter="\t")
np.savetxt(os.path.join(os.getcwd(),FolderName,FieldY_fromElectron), arElecEy1, newline="\n", delimiter="\t")
np.savetxt(os.path.join(os.getcwd(),FolderName,FieldX_fromBeam), arBeamEx1, newline="\n", delimiter="\t")
np.savetxt(os.path.join(os.getcwd(),FolderName,FieldY_fromBeam), arBeamEy1, newline="\n", delimiter="\t")
print('done in', round(time.time() - t0), 's')

#***********Plotting the Calculation Results
unitsIntPlot = ['m', 'm', 'ph/s/.1%bw/mm^2']
#uti_plot2d1d(arI0, [mesh0.xStart, mesh0.xFin, mesh0.nx], [mesh0.yStart, mesh0.yFin, mesh0.ny], labels=('Horizontal position', 'Vertical position', 'Intensity Before Lens'), units=unitsIntPlot)
uti_plot2d1d(arElecEx, [mesh0.xStart, mesh0.xFin, mesh0.nx], [mesh0.yStart, mesh0.yFin, mesh0.ny], labels=('Horizontal position', 'Vertical position', 'One Electron, Hor. polar.'), units=unitsIntPlot)
uti_plot2d1d(arElecEy, [mesh0.xStart, mesh0.xFin, mesh0.nx], [mesh0.yStart, mesh0.yFin, mesh0.ny], labels=('Horizontal position', 'Vertical position', 'One Electron, Vrt. polar.'), units=unitsIntPlot)
uti_plot2d1d(arBeamEx, [mesh0.xStart, mesh0.xFin, mesh0.nx], [mesh0.yStart, mesh0.yFin, mesh0.ny], labels=('Horizontal position', 'Vertical position', 'Whole Beam, Hor. polar.'), units=unitsIntPlot)
uti_plot2d1d(arBeamEy, [mesh0.xStart, mesh0.xFin, mesh0.nx], [mesh0.yStart, mesh0.yFin, mesh0.ny], labels=('Horizontal position', 'Vertical position', 'Whole Beam, Vrt. polar.'), units=unitsIntPlot)
uti_plot_show() #show all graphs (blocks script execution; close all graph windows to proceed)
