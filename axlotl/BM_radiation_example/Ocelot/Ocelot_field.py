from bmrad import *
import matplotlib.pyplot as plt
import numpy as np
import time
import csv
import multiprocessing
import os

#***********Data Folder and File Names
strFolderName = 'TwoD_Interferometer' #example data sub-folder name
strFieldFile_x = 'Ocelot_Field_X.txt' #file name for horizontal SR data
strFieldFile_y = 'Ocelot_Field_Y.txt' #  file name for vertical SR data
outfile_x = os.path.join(os.path.dirname(os.path.realpath(__file__)), strFolderName, strFieldFile_x)
outfile_y = os.path.join(os.path.dirname(os.path.realpath(__file__)), strFolderName, strFieldFile_y)


#queue = multiprocessing.Queue()
output_lock = multiprocessing.Lock()
worker_count = multiprocessing.cpu_count()

#defining a "worker" for further paralellization of the calculation
def worker(num, xr, yr, xo, yo, xc, yc, xpc, ypc, sgx, sgxp, sgy, sgyp, Ex, Ey, Lx, Ly):
    output_lock.acquire()
    print('Proc. number', num, 'has been started')
    output_lock.release()
    for i in range(int(xr[0]), int(xr[1])):
        for j in range(int(yr[0]), int(yr[1])):
            Ex[i,j], Ey[i,j] = bm_e_a(xc, yc, xpc, ypc, xo[i], yo[j], 14.2, sgx, sgxp, sgy, sgyp, p_en=2.748, R=191.73)
    Lx[num]=Ex
    Ly[num]=Ey

#beam parameters
sigx = 169.e-6
sigy = 10.e-6
sigxp = 0.000125
sigyp = 1.e-8

#Grids through horizontal and vertical electron position and its deflection in beam
yc = np.linspace(-0.*sigy,0.*sigy,1)
xc = np.linspace(-0.*sigx,0.*sigx,1)
ypc = np.linspace(-0.*sigyp,0.*sigyp,1)
xpc = np.linspace(-0.*sigxp,0.*sigxp,1)

#defining horizontal and vertical grid for the finite field
Ny = 4*worker_count
Nx = 4*worker_count
yo = np.linspace(-0.1e-1,0.1e-1, Ny)
xo = np.linspace(-0.1e-1,0.1e-1, Nx)
#empty matrices for Ex and Ey field components
Ex = np.zeros([len(xo),len(yo)], dtype='complex')
Ey = np.zeros([len(xo),len(yo)], dtype='complex')

print('Calculating SR field ...', end='')
#Begin of paralellization process
if __name__ == '__main__':
    mngr = multiprocessing.Manager()
    Ex_full = mngr.list([[] for i in range(worker_count)])
    Ey_full = mngr.list([[] for i in range(worker_count)])
    t0 = time.time()

    N=(len(yo))/worker_count

    yrange=[[]]
    for i in range(worker_count):
        if i==0:
            yrange[i]=[0, N]
        else:
            yrange.append([i*N, (i+1)*N])
    
    xrange = [0, len(xo)]
    jobs = []

    #******Creating list of processes
    for i in range(worker_count):
        proc = multiprocessing.Process(target=worker,
                                       args=(i, xrange, yrange[i],
                                             xo, yo, xc, yc, xpc, ypc,
                                             sigx, sigxp, sigy, sigyp,
                                             Ex, Ey, Ex_full, Ey_full))
        jobs.append(proc)
        proc.start()

    #******Finishing all the processes
    for job in jobs:
        job.join()

    print(' is done for ', round(time.time()-t0), 's')
    
    Lx = np.abs(np.sum(Ex_full, axis=0).T)**2
    Ly = np.abs(np.sum(Ey_full, axis=0).T)**2

    #*****************************************Plotting horizontal polarization
    ax1 = figure(1, figsize=(15,6)).add_subplot(131)
    p1 = ax1.imshow((Lx), extent = [xo[0]*1e+3,xo[-1]*1e+3,yo[0]*1e+3,yo[-1]*1e+3], aspect = 1.4*xo[0]/yo[0])
    p1.set_cmap('gray')
    ax1.set_xlabel('X [mm]')
    ax1.set_ylabel('Y [mm]')

    ax2 = figure(1).add_subplot(132)
    p2 = ax2.plot(yo*1e+3, Lx.T[round(Nx/2)])
    ax2.set_xlabel('Y [mm]')
    ax2.set_ylabel('arb. units')
    #ax2.set_title('Vertical')
    #print(L[6], '\n', end='')

    ax3 = figure(1).add_subplot(133)
    p3 = ax3.plot(xo*1e+3, Lx[round(Ny/2)])
    ax3.set_xlabel('X, [mm]')
    ax3.set_ylabel('arb. units')
    #ax3.set_title('Horizontal')
    #print(L.T[6], '\n', end='')

    #*****************************************Plotting vertical polarization
    ax4 = figure(2, figsize=(15,6)).add_subplot(131)
    p4 = ax4.imshow((Ly), extent = [xo[0]*1e+3,xo[-1]*1e+3,yo[0]*1e+3,yo[-1]*1e+3], aspect = 1.4*xo[0]/yo[0])
    p4.set_cmap('gray')
    ax4.set_xlabel('X [mm]')
    ax4.set_ylabel('Y [mm]')

    ax5 = figure(2).add_subplot(132)
    p5 = ax5.plot(yo*1e+3, Ly.T[round(Nx/2)])
    ax5.set_xlabel('Y [mm]')
    ax5.set_ylabel('arb. units')
    #ax2.set_title('Vertical')
    #print(L[6], '\n', end='')

    ax6 = figure(2).add_subplot(133)
    p6 = ax6.plot(xo*1e+3, Ly[round(Ny/2)])
    ax6.set_xlabel('X, [mm]')
    ax6.set_ylabel('arb. units')
    #ax3.set_title('Horizontal')
    #print(L.T[6], '\n', end='')

    plt.show()
    #figure(1).add_subplot(111).imshow(np.abs(Field)).set_cmap('gray')
    #plt.show()

    #Saving Ex and Ey matrices
    print('Saving field to a file ...', end='')
    np.savetxt(outfile_x, np.sum(Ex_full, axis=0).T, newline="\n", delimiter="\t")
    np.savetxt(outfile_y, np.sum(Ey_full, axis=0).T, newline="\n", delimiter="\t")
    print(' done')
    
