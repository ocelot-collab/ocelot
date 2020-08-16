'''
interface to srw
'''

import os
import string
from array import *

import sys
import numpy as np
import copy

from ocelot.cpbd.optics import *
from ocelot.cpbd.beam import *
from ocelot.cpbd.elements import *
from ocelot.rad.undulator_params import UndulatorParameters

from lib.srwlib import *
#from lib.srwlib import *
import os 


class SRRunParameters:
    def __init__(self):

        self.polarization = "POLARIZATION_TOTAL"
        self.method = 1
        self.methods = {"MANUAL":0, "AUTO_UNDULATOR":1, "AUTO_WIGGLER":2}
        self.nruns  = "1"
        self.nTrajectoryPoints = 1000000
        self.ctEnd = 200.0
        self.precision = 1.e-5

        
# wrapper class for srwlwavefront with getters and setters
class WaveFront:
    def __init__(self, srwlwfr=None):
        if srwlwfr != None:
            self.wfr = srwlwfr
       
    def getEx(self):
        pass
    def getEy(self):
        pass

def readField(filePath, sCom):
    f = open(filePath, 'r')
    f.readline() #1st line: just pass
    xStart = float(f.readline().split(sCom, 2)[1]) #2nd line: initial X position [m]; it will not actually be used
    xStep = float(f.readline().split(sCom, 2)[1]) #3rd line: step vs X [m]
    xNp = int(f.readline().split(sCom, 2)[1]) #4th line: number of points vs X
    yStart = float(f.readline().split(sCom, 2)[1]) #5th line: initial Y position [m]; it will not actually be used
    yStep = float(f.readline().split(sCom, 2)[1]) #6th line: step vs Y [m]
    yNp = int(f.readline().split(sCom, 2)[1]) #7th line: number of points vs Y
    zStart = float(f.readline().split(sCom, 2)[1]) #8th line: initial Z position [m]; it will not actually be used
    zStep = float(f.readline().split(sCom, 2)[1]) #9th line: step vs Z [m]
    zNp = int(f.readline().split(sCom, 2)[1]) #10th line: number of points vs Z
    totNp = xNp*yNp*zNp
    locArBx = array('d', [0]*totNp)
    locArBy = array('d', [0]*totNp)
    locArBz = array('d', [0]*totNp)
    for i in range(totNp):
        curLineParts = f.readline().split('\t')
        locArBx[i] = float(curLineParts[0])
        locArBy[i] = float(curLineParts[1])
        locArBz[i] = float(curLineParts[2])
    f.close()
    xRange = xStep
    if xNp > 1: xRange = (xNp - 1)*xStep
    yRange = yStep
    if yNp > 1: yRange = (yNp - 1)*yStep
    zRange = zStep
    if zNp > 1: zRange = (zNp - 1)*zStep
    return SRWLMagFld3D(locArBx, locArBy, locArBz, xNp, yNp, zNp, xStep*(xNp - 1), yStep*(yNp - 1), zStep*(zNp - 1), 1)


def saveTrajectory(traj, filePath):
    f = open(filePath, 'w')
    f.write('#ct [m], X [m], BetaX [rad], Y [m], BetaY [rad], Z [m], BetaZ [m]\n')
    ctStep = 0
    if traj.np > 0:
        ctStep = (traj.ctEnd - traj.ctStart)/(traj.np - 1)
    ct = traj.ctStart
    for i in range(traj.np):
        f.write(str(ct) + '\t' + repr(traj.arX[i]) + '\t' + repr(traj.arXp[i]) + '\t' + repr(traj.arY[i]) + '\t' + repr(traj.arYp[i]) + '\t' + repr(traj.arZ[i]) + '\t' + repr(traj.arZp[i]) + '\n')        
        ct += ctStep
    f.close()

def readTrajectory(filePath):

    #traj = SRWLPrtTrj()
    
    traj = {}
    
    npTraj = 0

    f=open(filePath,'r')
    f.readline()
    for line in f:
        if len(line.strip()) > 1:
            npTraj += 1
  
    #traj.allocate(npTraj)

    traj['x'] = array('d', [0]*npTraj)
    traj['xp'] = array('d', [0]*npTraj)
    traj['y'] = array('d', [0]*npTraj)
    traj['yp'] = array('d', [0]*npTraj)
    traj['z'] = array('d', [0]*npTraj)
    traj['zp'] = array('d', [0]*npTraj)


    f.close()    
    f=open(filePath,'r')
    
    f.readline()
    i = 0
    for line in f:
        if i < npTraj:
            #ct, traj.arX[i], traj.arXp[i], traj.arY[i], traj.arYp[i], traj.arZ[i], traj.arZp[i] = map(float, line.strip().split('\t'))
            ct, traj['x'][i], traj['xp'][i], traj['y'][i], traj['yp'][i], traj['z'][i], traj['zp'][i] = map(float, line.strip().split('\t'))
            
            if i == 0:
                traj['ctSrart'] = ct
            if i == npTraj-1:
                traj['ctEnd'] = ct

            i += 1
 

    f.close()
    return traj


def parseParameter(token, parameter):
    if token.startswith(parameter+"="):
        return float(token.replace(parameter+"=",""))
    

def readIntensity(filePath, type='spectral'):
    
    #print type
    
    intens = Intensity()
    s = []
    val = []
    
    intens.nx=1
    intens.ny=1
    
    intens.xmin = 0
    intens.xmax = 0 
    intens.ymin = 0
    intens.ymax = 0
    
    w0=0.0
    wend = 1.0
    nPointsSpectrum = 10
    dw=(wend - w0)/nPointsSpectrum
    
    startData = True
    
    idx = 0
    
    f=open(filePath,'r')
    for line in f: 
        if line.strip().startswith("#"):
            pass
        elif line.strip().startswith("@"):
          
            token = line.strip().split()[0]
            #print token
          
            if token.startswith("@startPhotonEnergy="):
                w0=parseParameter(token, "@startPhotonEnergy")
            if token.startswith("@endPhotonEnergy="):
                wend=parseParameter(token, "@endPhotonEnergy")
            if token.startswith("@photonEnergyPoints="):
                nPointsSpectrum=int(parseParameter(token, "@photonEnergyPoints"))
            if token.startswith("@xPoints="):
                intens.nx=parseParameter(token, "@xPoints")
            if token.startswith("@yPoints="):
                intens.ny=parseParameter(token, "@yPoints")
            if token.startswith("@startX="):
                intens.xmin=parseParameter(token, "@startX")
            if token.startswith("@endX="):
                intens.xmax=parseParameter(token, "@endX")
            if token.startswith("@startY="):
                intens.ymin=parseParameter(token, "@startY")
            if token.startswith("@endY="):
                intens.ymax=parseParameter(token, "@endY")
             
        else:
            if startData:
                if nPointsSpectrum > 1:
                    dw=(wend - w0)/(nPointsSpectrum-1)
                w = w0
                startData = False

                for i in range(0,nPointsSpectrum):
                    w = round(w0 + i*dw,4)
                    intens.intensity[w] = np.zeros([intens.nx, intens.ny])

                w = w0

                #print 'spectral points', intens.intensity

                iSpec = 0

            try:
                if type == 'spectrum':
                    val.append(float(line.strip()))
                    s.append(w)
                    w = w + dw

                if type == 'intensity':

                    j = int(idx % intens.nx)
                    i = int(idx / intens.nx)

                    intens.intensity[round(w0+iSpec*dw,4)][i,j] = float(line.strip())

                    if iSpec == nPointsSpectrum-1:
                        idx += 1

                    iSpec = (iSpec + 1) % nPointsSpectrum
                  

            except ValueError:
                pass
    
    intens.spectrum[(0,0)] = [s,val]
    
    return intens



# operations on intenisties
def diff(int1,int2):
    int = Intensity(int1)
    int.intensitySlice = int1.intensitySlice - int2.intensitySlice 
    print (int.intensitySlice)
    return int

def sum(int1,int2):
    int = Intensity()
    int.intensitySlice = int1.intensitySlice + int2.intensitySlice
    return int

def mean(int):
    return np.mean(int.intensitySlice)


    
def calculateSR_py(lat, beam, screen, runParameters):
    
    und = UndulatorParameters()

    e = lat.sequence[0] 
        
    #print e
    #print 'starting with', e.id, e.type
        
    pos = -e.l / 2.0
        
    if e.type == 'undulator':
        nw = int (e.nperiods)
        lw = float (e.lperiod)
        z0 = -0.5* lw*(nw + 4) # start of integration
        #print 'undulator z0=', z0

    else:
        l = float (e.l)
        z0 = -0.5* l 
        #print 'undulator z0=', z0



    fields = []
    arzc = []
    arxc = []
    aryc = []
    z = 0.0

    for e in lat.sequence:
        
        if e.type == "undulator":
            
            und.K = float(e.Kx)
            und.E = beam.E
            und.lw = float(e.lperiod) 
            und.Nw = float(e.nperiods)    
            und.recalculate()
    
            numPer = e.nperiods
            undPer = e.lperiod #Period Length [m]
            Bx = 0.0 #Peak Horizontal field [T]
            By = und.B #Peak Vertical field [T]
            phBx = 0 #Initial Phase of the Horizontal field component
            phBy = 0 #Initial Phase of the Vertical field component
            sBx = -1 #Symmetry of the Horizontal field component vs Longitudinal position
            sBy = 1 #Symmetry of the Vertical field component vs Longitudinal position

            fields.append(SRWLMagFldU([SRWLMagFldH(1, 'v', By, phBy, sBy, 1), 
                                       SRWLMagFldH(1, 'h', Bx, phBx, sBx, 1)], undPer, numPer))
            
            
            arzc.append(pos + e.l / 2.0)
            arxc.append(e.dx)
            aryc.append(e.dy)
            
        if e.type == "quadrupole":
            
            kn = beam.E * e.k1 / 0.2998
                        
            fields.append(SRWLMagFldM(kn, 2, 'n', e.l))
            arzc.append(pos + e.l / 2.0)
            arxc.append(e.dx)
            aryc.append(e.dy)

        elif e.type == 'rbend' or e.type == 'sbend' or e.type == 'bend':
                        
            Bx = 0.0
            By = -3.3356 * e.angle * beam.E / e.l 
            Bz = 0.0
            
            fields.append(SRWLMagFldM(By, 1, 'n', e.l))
            arzc.append(pos + e.l / 2.0)
            arxc.append(e.dx)
            aryc.append(e.dy)

        
        pos += e.l

    magFldCnt = SRWLMagFldC(fields, array('d',arxc), array('d',aryc), array('d',arzc)) #Container of all Field Elements

    zstart = z0
    zend = runParameters.ctEnd 
    #print 'z:', zstart, zend
    #print screen.nx, screen.ny

    #***********Precision
    meth = 1 #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
    relPrec = runParameters.precision #relative precision
    zStartInteg = zstart #longitudinal position to start integration (effective if < zEndInteg)
    zEndInteg = zend #longitudinal position to finish integration (effective if > zStartInteg)
    npTraj = runParameters.nTrajectoryPoints
    sampFactNxNyForProp = 0 #sampling factor for adjusting nx, ny (effective if > 0)
    arPrecPar = [meth, relPrec, zStartInteg, zEndInteg, npTraj, 0, sampFactNxNyForProp]


    #**********************Trajectory for the reference particle
    
    part = SRWLParticle()
    part.x = beam.x 
    part.y = beam.y
    part.xp = beam.xp
    part.yp = beam.yp
    part.gamma = beam.E/0.51099890221e-03 #Relative Energy
    part.relE0 = 1 #Electron Rest Mass
    part.nq = -1 #Electron Charge
    part.z = z0
    
    partTraj = SRWLPrtTrj()
    partTraj.partInitCond = part
    partTraj.allocate(npTraj)
    partTraj.ctStart = 0 #Start Time for the calculation
    partTraj.ctEnd = runParameters.ctEnd

    print ('Performing trajectory calculation ... ')
    partTraj = srwl.CalcPartTraj(partTraj, magFldCnt, arPrecPar)
    print ('done')


    #**********************Field calculation

    elecBeam = SRWLPartBeam()
    elecBeam.Iavg = beam.I #Average Current [A]
    elecBeam.partStatMom1.x = beam.x #Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
    elecBeam.partStatMom1.y = beam.y
    elecBeam.partStatMom1.z = zstart #Initial Longitudinal Coordinate (set before the ID)
    elecBeam.partStatMom1.xp = beam.xp #Initial Relative Transverse Velocities
    elecBeam.partStatMom1.yp = beam.yp
    elecBeam.partStatMom1.gamma = beam.E/0.51099890221e-03 #Relative Energy
    
    wfr1 = SRWLWfr() #For spectrum vs photon energy

    wfr1.allocate(screen.num_energy, screen.nx, screen.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
    wfr1.mesh.zStart = screen.z #Longitudinal Position [m] at which SR has to be calculated
    wfr1.mesh.eStart = screen.start_energy #Initial Photon Energy [eV]
    wfr1.mesh.eFin = screen.end_energy #Final Photon Energy [eV]
    wfr1.mesh.xStart = -screen.size_x / 2.0 + screen.x #Initial Horizontal Position [m]
    wfr1.mesh.xFin = screen.size_x / 2.0 + screen.x #Final Horizontal Position [m]
    wfr1.mesh.yStart = -screen.size_y / 2.0 + screen.y#Initial Vertical Position [m]
    wfr1.mesh.yFin = screen.size_y / 2.0 + screen.y #Final Vertical Position [m]
    wfr1.partBeam = elecBeam

    
    print ('Performing Electric Field (spectrum vs photon energy) calculation ... ')
    arPrecPar = [runParameters.method, relPrec, zStartInteg, zEndInteg, npTraj, 0, sampFactNxNyForProp]
    srwl.CalcElecFieldSR(wfr1, 0, magFldCnt, arPrecPar)
    print ('done')
    print ('   Extracting Intensity from calculated Electric Field ... ')
    arI1 = array('f', [0]*wfr1.mesh.ne*wfr1.mesh.nx*wfr1.mesh.ny)
    srwl.CalcIntFromElecField(arI1, wfr1, 6, 0, 6, wfr1.mesh.eStart, wfr1.mesh.xStart, wfr1.mesh.yStart)
    print ('done')
    
    intensity = np.array(arI1).reshape(screen.nx, screen.ny,screen.num_energy)
    
    return (partTraj, intensity)
    