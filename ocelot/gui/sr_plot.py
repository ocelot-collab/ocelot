__author__ = 'Sergey Tomin'

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
from matplotlib import cm


def show_flux(screen, show='Total', xlim=(0,0), ylim =(0,0), file_name=None, unit="mm", title=None, nfig=1, grid=True,
              xlog=False, ylog=False):
    """
    Function to plot radiation flux density 1D and 3D

    :param screen: Screen class
    :param show: str, show polarization component, can be  "Total", "Sigma" or "Pi"
    :param xlim: (0,0),  xlimits, ignored if (0,0)
    :param ylim: (0,0),  ylimits, ignored if (0,0)
    :param file_name: str, if None flux will be shown otherwise save to file
    :param unit: str, units for spatial axis, can be "mm" or "mrad",
    :param title: figure title
    :param nfig: figure number
    :param grid: True, if False not show the grid
    :param xlog: False, log scale for x-axis
    :param ylog: False, log scale for y-axis
    :return:
    """
    if show == 'Total':
        data = screen.Total
    elif show == 'Sigma':
        data = screen.Sigma
    else:
        data = screen.Pi

    if screen.nx == 1 or screen.ny == 1:
        if screen.nx == 1 and screen.ny == 1:
            X = screen.Eph
            xlabel = r'$E_{ph}$, $eV$'
            status = "spectrum"
        elif screen.nx == 1:
            X = screen.Yph
            xlabel = r'$Y$, $mm$'
            if unit == "mrad":
                xlabel = r'$Y$, $mrad$'
            status = "spatial"
        else:
            X = screen.Xph
            xlabel = r'$X$, $mm$'
            if unit == "mrad":
                xlabel = r'$X$, $mrad$'
            status = "spatial"

        D1(data, X, distance=screen.Distance, xlabel=xlabel, xlim=xlim, ylim=ylim,  file_name=file_name,
           unit=unit, status=status, title=title, nfig=nfig, grid=grid, xlog=xlog, ylog=ylog)
    else:
        if screen.ne != 1:
            print (" ******** Can not display 4D plot. Change number of points Screen.num_energy = 1 ! *********** ")
            return
        D3(screen, data, distance=screen.Distance, file_name=file_name, unit=unit, title=title, nfig=nfig)


def D1(data, X, distance, xlabel, xlim, ylim,  file_name, unit, status, xlog=False, ylog=False, title=None, nfig=1, grid=True):
    # distance in [mm]
    if unit == "mrad":
        data = data*distance*distance*1e-6
    if unit == "mrad" and status == "spatial":
        X = X/distance*1e3

    maxS = max(data)
    index = np.where(data== max(data))[0][0]
    energy = X[index]

    fig = plt.figure(nfig)

    if title is not None:
        fig.suptitle(title, fontsize=15)

    ax = fig.add_subplot(111)

    ax.plot(X, data)

    if xlim != (0, 0):
        ax.set_xlim(xlim)
    if ylim != (0, 0):
        ax.set_ylim(ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$I$, $\frac{ph}{sec \cdot mm^2 10^{-3}BW}$")
    if unit == "mrad":
        ax.set_ylabel(r"$I$, $\frac{ph}{sec \cdot mrad^2 10^{-3}BW}$")
    ax.grid(grid)
    if ylog: ax.set_yscale('log')
    if xlog: ax.set_xscale('log')
    # ax.annotate('$\epsilon_1 = ' + str(int(energy*10)/10.) +'$', xy=(0.9, 0.85),
    #            xycoords='axes fraction',
    #            horizontalalignment='right', verticalalignment='top',
    #            fontsize=20)

    power = np.floor(np.log10(maxS))
    intensity = np.around(maxS*10**(-power), 2)*10**power
    ax.annotate('I = ' + str("{0:.2E}".format(intensity)) , xy=(0.9, 0.93),
               xycoords='axes fraction',
               horizontalalignment='right', verticalalignment='top',
               fontsize=15)

    if file_name != None:
        figg = plt.gcf()
        k_size = 1.4
        figg.set_size_inches( (4*k_size, 3.01*k_size) )
        figg.savefig(file_name)
    else:
        plt.show()


def D3(screen,Data, distance, file_name = None , unit = "mm", title=None, nfig=1):
    #print " showme.any = ", np.shape(Data)
    X,Y = np.meshgrid(screen.Xph, screen.Yph)
    if unit == "mrad":
        Data = Data*distance*distance*1e-6
        X = X/distance*1e6
        Y = Y/distance*1e6
    fig = plt.figure(nfig)
    if title is not None:
        plt.title(title)
    ax = fig.add_subplot(111, projection='3d')


    #print " showme.any = ", np.shape(X)
    #print " showme.any = ", np.shape(Y)
    data = np.zeros((screen.ny, screen.nx))
    for j in range(screen.ny):
        for i in range(screen.nx):
            data[j,i] = Data[screen.nx*j + i]
    ax.plot_surface(X, Y, data, rstride=1, cstride=1, cmap=cm.jet)
    #ax.set_zlim3d(0, 1)

    if unit == "mrad":
        ax.set_xlabel(r'$\theta_x$, $\mu rad$')
        ax.set_ylabel(r'$\theta_y$, $\mu rad$')
        ax.set_zlabel(r"$I$, $\frac{ph}{s \cdot mrad^2 10^{-3}BW}$")
    else:
        ax.set_xlabel(r'$X$, $mm$')
        ax.set_ylabel(r'$Y$, $mm$')
        ax.set_zlabel(r"$I$, $\frac{ph}{s\cdot mm^2 10^{-3}BW}$")

    #ax.set_xticks([])
    if file_name != None:
        figg = plt.gcf()
        k_size = 1.7
        figg.set_size_inches( (4*k_size, 3.01*k_size) )
        figg.savefig(file_name)
    else:
        plt.show()
    #plt.show()

def plot3D_data(data, x = None, y = None):
    if x != None and y != None:
        X,Y = np.meshgrid(x,y)
    else:
        print( np.shape(data))
        X,Y = np.meshgrid(np.arange(np.shape(data)[1]), np.arange(np.shape(data)[0]))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, data, rstride=1, cstride=1, cmap=cm.jet)
    plt.show()