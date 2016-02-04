'''
user interface for viewing genesis simulation results
'''

import sys, os, csv

 
import matplotlib
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
# from matplotlib.collections import PatchCollection
# import matplotlib.patches as mpatches
# import matplotlib.path as mpath
import matplotlib.pyplot as plt
# from ocelot.cpbd.optics import *
import numpy as np
from numpy import matlib
# from ocelot.common.math_op import *

from pylab import * #tmp
params = {'backend': 'ps', 'axes.labelsize': 15, 'font.size': 16, 'legend.fontsize': 24, 'xtick.labelsize': 19,  'ytick.labelsize': 19, 'text.usetex': True}
# params = {'backend': 'ps', 'axes.labelsize': 18, 'text.fontsize': 16, 'legend.fontsize': 24, 'xtick.labelsize': 32,  'ytick.labelsize': 32, 'text.usetex': True}
rcParams.update(params)
rc('text', usetex=True) # required to have greek fonts on redhat


h = 4.135667516e-15
c = 299792458.0

max_yticks = 7

def gen_outplot_e(g, figsize=(8,10), legend = True, fig_name = None):

    import matplotlib.ticker as ticker

    font_size = 1
    if fig_name is None:
        if g.filename is '':
            fig = plt.figure('Electrons')
        else:
            fig = plt.figure('Electrons '+g.filename)
    else:
        fig = plt.figure(fig_name)
    
    fig.set_size_inches(figsize,forward=True)
    plt.rc('axes', grid=True)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    # left, width = 0.1, 0.85
    plt.clf()
    #
    # rect1 = [left, 0.75, width, 0.20]
    # rect2 = [left, 0.55, width, 0.20]
    # rect3 = [left, 0.35, width, 0.20]
    # rect4 = [left, 0.15, width, 0.20]

    #
    ax_und=fig.add_subplot(4, 1, 1)
    ax_und.clear()
    ax_sizepos=fig.add_subplot(4, 1, 2,sharex=ax_und)
    ax_sizepos.clear()
    ax_energy=fig.add_subplot(4, 1, 3,sharex=ax_und)
    ax_energy.clear()
    ax_bunching=fig.add_subplot(4, 1, 4,sharex=ax_und)
    ax_bunching.clear()


    # for ax in ax_sizepos, ax_energy, ax_und, ax_bunching:

    # ax_sizepos.grid(True)
    # ax_und.grid(True)
    # ax_energy.set_yticks([])
    # ax_energy.grid(True)


    # ax_und = fig.add_axes(rect1)
    # ax_sizepos = fig.add_axes(rect2, sharex=ax_und)  #left, bottom, width, height
    # ax_energy = fig.add_axes(rect3, sharex=ax_und)
    # ax_bunching = fig.add_axes(rect4, sharex=ax_und)
    for ax in ax_sizepos, ax_energy, ax_und, ax_bunching:
        if ax!=ax_bunching:
            for label in ax.get_xticklabels():
                label.set_visible(False)


    # for tick in ax.yaxis.get_major_ticks():
    #     tick.label.set_fontsize(14)
    #     # specify integer or one of preset strings, e.g.
    #     #tick.label.set_fontsize('x-small')
    #     tick.label.set_rotation('vertical')




    #
    fig.subplots_adjust(hspace=0)
    # # beta_x = [p.beta_x for p in tws] # list(map(lambda p:p.beta_x, tws))
    # # beta_y = [p.beta_y for p in tws] #list(map(lambda p:p.beta_y, tws))
    # # S = [p.s for p in tws] #list(map(lambda p:p.s, tws))
    # #plt.plot(S, beta_x)

#    t_array = np.linspace(g.t[0], g.t[-1], len(g.t))
#    s_array = t_array*c*1.0e-15
#    ncar = g('ncar')
#    dgrid = g('dgrid')
    
#    xlamds = float(g('xlamds'))
#    nslice = len(g.spec)
#    zsep  = float(g('zsep'))
#    zstop  = float(g('zstop'))
#    
#    srange = nslice*zsep*xlamds     #range in s
#    dk = 2*np.pi/srange
#    ds = zsep*xlamds
#    
#    #w_l_m  =  g('xlamds')
#    w_l_ev = h * c / g('xlamds')
    
    # max_pw_list = np.zeros(run_end+1)
    # max_argpw_list = np.zeros(run_end+1)
    # max_sp_list = np.zeros(run_end+1)
    # max_argsp_list = np.zeros(run_end+1)
    
    #
    
    # plot_lattice(g, font_size)
    
    # plot_betas(ax_sizepos, S, beta_x, beta_y, font_size)
    
    # plot_elems(ax_energy, lat, s_point = S[0], legend = legend, y_scale=0.8) # plot elements

 

    #sys.exit()

    # ax_und
    # ax_sizepos
    # ax_energy
    # ax_spread

    ax_und.plot(g.z, g.aw, 'b-',linewidth=1.5)
    ax_und.set_ylabel('K')

    ax_quad = ax_und.twinx()
    ax_quad.plot(g.z, g.qfld, 'r-',linewidth=1.5)
    ax_quad.set_ylabel('Quad')

    # print len(g.z),len(g.xrms[0,:]),len(np.mean(g.yrms,axis=0))

    #sys.exit()
    ax_sizepos.plot(g.z, np.mean(g.xrms,axis=0)*1e6, 'g-',g.z, np.mean(g.yrms,axis=0)*1e6, 'b-')
    ax_sizepos.set_ylabel('$\sigma_{x,y}$ [$\mu$m]')


    ax_energy.plot(g.z, np.average(g.el_energy*0.511e-3, weights=g.I, axis=0), 'b-',linewidth=1.5)
    ax_energy.set_ylabel('E [GeV]')
    ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    #ax_energy.plot(g.z, np.mean(g.energy_GeV,axis=0), 'g-')
    ax_spread = ax_energy.twinx()
    ax_spread.plot(g.z, np.average(g.el_e_spread*0.511e-3*1000, weights=g.I, axis=0), 'r--', g.z, np.amax(g.el_e_spread*0.511e-3*1000, axis=0), 'r-')
    ax_spread.set_ylabel('$\sigma_E$ [MeV]')

    number_ticks=4


    #ax_bunching.plot(g.z, np.average(g.bunching, weights=g.I[1:], axis=0), 'k--', g.z, np.amax(g.bunching, axis=0), 'r-')    
    ax_bunching.plot(g.z, np.average(g.bunching, weights=g.I, axis=0), 'k--', g.z, np.amax(g.bunching, axis=0), 'g-')
    ax_bunching.set_ylabel('Bunching')

    ax_bunching.set_xlabel('s [$\mu$m]')
    x1,x2,y1,y2 = ax_sizepos.axis()
    ax_sizepos.axis([x1,x2,0,y2])
    x1,x2,y1,y2 = ax_spread.axis()
    ax_spread.axis([x1,x2,0,y2])

    ax_und.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_quad.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_sizepos.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_spread.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_bunching.yaxis.major.locator.set_params(nbins=number_ticks)

    # yloc = plt.MaxNLocator(max_yticks)
    # # print yloc
    # ax_sizepos.yaxis.set_major_locator(yloc)
    # ax_energy.yaxis.set_major_locator(yloc)
    # ax_und.yaxis.set_major_locator(yloc)
    # ax_bunching.yaxis.set_major_locator(yloc)

    # ax_energy.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1e'))

    plt.xlim(g.z[0], g.z[-1])

    fig.subplots_adjust(top=0.95, bottom=0.1, right=0.85, left=0.15)

    ax_und.tick_params(axis='y', which='both', colors='b')
    ax_und.yaxis.label.set_color('b')    
    ax_quad.tick_params(axis='y', which='both', colors='r')
    ax_quad.yaxis.label.set_color('r') 
    ax_energy.tick_params(axis='y', which='both', colors='b')
    ax_energy.yaxis.label.set_color('b')    
    ax_spread.tick_params(axis='y', which='both', colors='r')
    ax_spread.yaxis.label.set_color('r') 



    plt.show()


def gen_outplot_ph(g, figsize=(8, 10), legend = True, fig_name = None):
#    max_yticks = 7
    import matplotlib.ticker as ticker
    
    
    font_size = 1
    if fig_name is None:
        if g.filename is '':
            fig = plt.figure('Photons')
        else:
            fig = plt.figure('Photons '+g.filename)
    else:
        fig = plt.figure(fig_name)

    fig.set_size_inches(figsize,forward=True)

    plt.rc('axes', grid=True)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    # left, width = 0.1, 0.85
    plt.clf()


    #
    # rect1 = [left, 0.75, width, 0.20]
    # rect2 = [left, 0.55, width, 0.20]
    # rect3 = [left, 0.35, width, 0.20]
    # rect4 = [left, 0.15, width, 0.20]

    #
    ax_pow=fig.add_subplot(3, 1, 1)
    ax_pow.clear()
    ax_spectrum=fig.add_subplot(3, 1, 2,sharex=ax_pow)
    ax_spectrum.clear()
    ax_size=fig.add_subplot(3, 1, 3,sharex=ax_pow)
    ax_size.clear()


    # for ax in ax_sizepos, ax_energy, ax_und, ax_bunching:

    # ax_sizepos.grid(True)
    # ax_und.grid(True)
    # ax_energy.set_yticks([])
    # ax_energy.grid(True)


    # ax_und = fig.add_axes(rect1)
    # ax_sizepos = fig.add_axes(rect2, sharex=ax_und)  #left, bottom, width, height
    # ax_energy = fig.add_axes(rect3, sharex=ax_und)
    # ax_bunching = fig.add_axes(rect4, sharex=ax_und)
    for ax in ax_pow, ax_spectrum, ax_size:
        if ax!=ax_size:
            for label in ax.get_xticklabels():
                label.set_visible(False)


    # for tick in ax.yaxis.get_major_ticks():
    #     tick.label.set_fontsize(14)
    #     # specify integer or one of preset strings, e.g.
    #     #tick.label.set_fontsize('x-small')
    #     tick.label.set_rotation('vertical')




    #
    fig.subplots_adjust(hspace=0)
    # # beta_x = [p.beta_x for p in tws] # list(map(lambda p:p.beta_x, tws))
    # # beta_y = [p.beta_y for p in tws] #list(map(lambda p:p.beta_y, tws))
    # # S = [p.s for p in tws] #list(map(lambda p:p.s, tws))
    # #plt.plot(S, beta_x)

#    t_array = np.linspace(g.t[0], g.t[-1], len(g.t))
#    s_array = t_array*c*1.0e-15
#    ncar = g('ncar')
#    dgrid = g('dgrid')
#
#    xlamds = float(g('xlamds'))
#    nslice = len(g.spec)
#    zsep  = float(g('zsep'))
#    zstop  = float(g('zstop'))
#
#    srange = nslice*zsep*xlamds     #range in s
#    dk = 2*np.pi/srange
#    ds = zsep*xlamds
#
#    #w_l_m  =  g('xlamds')
#    w_l_ev = h * c / g('xlamds')

    # max_pw_list = np.zeros(run_end+1)
    # max_argpw_list = np.zeros(run_end+1)
    # max_sp_list = np.zeros(run_end+1)
    # max_argsp_list = np.zeros(run_end+1)

    #

    # plot_lattice(g, font_size)

    # plot_betas(ax_sizepos, S, beta_x, beta_y, font_size)

    # plot_elems(ax_energy, lat, s_point = S[0], legend = legend, y_scale=0.8) # plot elements

#    print len(t_array)
#    print len(g.xrms)
#    print len(g.z)

    #sys.exit()

    # ax_und
    # ax_sizepos
    # ax_energy
    # ax_spread

    # ax_pow
    # ax_spectrum
    # ax_size

    ax_pow.plot(g.z, np.amax(g.p_int, axis=0), 'g-',linewidth=1.5)
    ax_pow.set_ylabel('P [W]')
    ax_pow.set_yscale('log')

    # outp.power.mean_S*inp.xlamds*inp.zsep*inp.nslice/c
    ax_en = ax_pow.twinx()
    ax_en.plot(g.z, np.mean(g.p_int,axis=0)*g('xlamds')*g('zsep')*g.nSlices/c, 'r--')
    ax_en.set_ylabel('E [J]')
    ax_en.set_yscale('log')

    n_pad=1
    # print len(g.z),len(g.xrms[0,:]),len(np.mean(g.yrms,axis=0))
    power=np.pad(g.p_mid, [(int(g.nSlices/2)*n_pad, (g.nSlices-(int(g.nSlices/2))))*n_pad, (0, 0)], mode='constant')
    phase=np.pad(g.phi_mid, [(int(g.nSlices/2)*n_pad, (g.nSlices-(int(g.nSlices/2))))*n_pad, (0, 0)], mode='constant')
    # spectrum = np.power(abs(fft(np.sqrt( np.array(g.p_mid)) * np.exp( 1.j* np.array(g.phi_mid) ) , axis=0)),2)#/sqrt(g.nSlices)
    spectrum = np.power(abs(fft(np.sqrt( np.array(power)) * np.exp( 1.j* np.array(phase) ) , axis=0)),2)/sqrt(g.nSlices)/np.power((2*g.leng/g('ncar')),2)/1e10
    e_0=1239.8/g('xlamds')/1e9
    # print e_0

    g.freq_ev1 = h * fftfreq(len(spectrum), d=g('zsep') * g('xlamds') / c)+e_0
    lamdscale=1239.8/g.freq_ev1
    #sys.exit()
#    print 'shape',spectrum.shape
#    spectrum_lamdpos=sum(np.matlib.repmat(lamdscale,spectrum.shape[1],1)*transpose(spectrum),axis=1)/np.sum(spectrum,axis=0)


    ax_spectrum.plot(g.z, np.amax(spectrum,axis=0), 'g-')
    ax_spectrum.set_ylabel('P_{max}($\lambda$)')
    ax_spectrum.set_yscale('log')

    #fix!!!
    # spectrum_lamdpos=sum(np.matlib.repmat(lamdscale,spectrum.shape[1],1)*transpose(spectrum),axis=1)/np.sum(spectrum,axis=0)
    # spectrum_width=sqrt(sum((np.power(np.matlib.repmat(lamdscale,spectrum.shape[1],1)-np.matlib.repmat(spectrum_lamdpos,spectrum.shape[1],1),2)*spectrum),1)/sum(spectrum,axis=0));
    # ax_spec_bandw = ax_spectrum.twinx()
    # ax_spec_bandw.plot(g.z, spectrum_lamdpos, 'g-')
    # fix and include!!!



    # print 'spec1', len(spectrum)
    # print 'freq_ev1', len(g.freq_ev1)


    #
    # fig2=plt.plot(fftshift(g.freq_ev1),fftshift(spectrum[:,-1],axes=0))
    # plt.show()
    # fig3=plt.plot(fftshift(lamdscale),fftshift(spectrum[:,-1],axes=0))
    # plt.show()

#    print g.r_size.shape
#    print g.p_int.shape

    g.p_int=np.amax(g.p_int)/1e6+g.p_int # nasty fix from division by zero
    ax_size.plot(g.z, np.average(g.r_size*2*1e6, weights=g.p_int, axis=0), 'b-')
    ax_size.set_ylabel('transverse [$\mu$m]')
    #ax_energy.plot(g.z, np.mean(g.energy_GeV,axis=0), 'g-')
    # ax_spread = ax_energy.twinx()
    # ax_spread.plot(g.z, np.average(g.e_spread_GeV, weights=g.I[1:], axis=0), 'b--', g.z, np.amax(g.e_spread_GeV, axis=0), 'b-')
    # ax_spread.set_ylabel('Energy spread ($\Delta$E)')

    number_ticks=4



    # ax_bunching.plot(g.z, np.average(g.bunching, weights=g.I[1:], axis=0), 'k--', g.z, np.amax(g.bunching, axis=0), 'r-')
    # ax_bunching.set_ylabel('Bunching')

    # ax_bunching.set_xlabel('s [m]')
    # x1,x2,y1,y2 = ax_sizepos.axis()
    # ax_sizepos.axis([x1,x2,0,y2])
    # x1,x2,y1,y2 = ax_spread.axis()
    # ax_spread.axis([x1,x2,0,y2])

# ax_pow
    # ax_spectrum
    # ax_size


    # ax_pow.yaxis.major.locator.set_params(nbins=number_ticks)
    # ax_en.yaxis.major.locator.set_params(nbins=number_ticks)
    # ax_spectrum.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_size.yaxis.major.locator.set_params(nbins=number_ticks)
    # ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)
    # ax_bunching.yaxis.major.locator.set_params(nbins=number_ticks)

    # yloc = plt.MaxNLocator(max_yticks)
    # # print yloc
    # ax_sizepos.yaxis.set_major_locator(yloc)
    # ax_energy.yaxis.set_major_locator(yloc)
    # ax_und.yaxis.set_major_locator(yloc)
    # ax_bunching.yaxis.set_major_locator(yloc)

    # ax_energy.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1e'))

    plt.xlim(g.z[0], g.z[-1])

    fig.subplots_adjust(top=0.95, bottom=0.1, right=0.85, left=0.15)


    ax_pow.tick_params(axis='y', which='both', colors='g')
    ax_pow.yaxis.label.set_color('g')    
    ax_en.tick_params(axis='y', which='both', colors='r')
    ax_en.yaxis.label.set_color('r') 
    ax_spectrum.tick_params(axis='y', which='both', colors='b')
    ax_spectrum.yaxis.label.set_color('b')    
    ax_size.tick_params(axis='y', which='both', colors='b')
    ax_size.yaxis.label.set_color('b') 


    plt.show()


def gen_outplot_z(g, figsize=(8, 10), legend = True, fig_name = None, z=inf):
#    max_yticks = 7
    import matplotlib.ticker as ticker

    if z==inf:
        print 'Showing profile parameters at the end of undulator'
        z=np.amax(g.z)

    if z>np.amax(g.z):
        print 'Z parameter too large, setting to the undulator end'
        z=np.amax(g.z)
    
    if z<np.amin(g.z):
        print 'Z parameter too small, setting to the undulator entrance'    
        z=np.amin(g.z)
        
    zi=np.where(g.z>=15)[0][0]
    z=g.z[zi];




    font_size = 1
    if fig_name is None:
        if g.filename is '':
            fig = plt.figure('Bunch profile')
        else:
            fig = plt.figure('Bunch profile '+g.filename)
    else:
        fig = plt.figure(fig_name)
    fig.set_size_inches(figsize,forward=True)
    plt.rc('axes', grid=True)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    # left, width = 0.1, 0.85
    
    plt.clf()
    

    #
    # rect1 = [left, 0.75, width, 0.20]
    # rect2 = [left, 0.55, width, 0.20]
    # rect3 = [left, 0.35, width, 0.20]
    # rect4 = [left, 0.15, width, 0.20]

    #
    ax_curr=fig.add_subplot(4, 1, 1)
    ax_curr.clear()
    ax_energy=fig.add_subplot(4, 1, 2,sharex=ax_curr)
    ax_energy.clear()    
    ax_phase=fig.add_subplot(4, 1, 3,sharex=ax_curr)
    ax_phase.clear()
    ax_spectrum=fig.add_subplot(4, 1, 4)
    ax_spectrum.clear()

    for ax in ax_curr, ax_phase, ax_spectrum, ax_energy:
        if ax!=ax_spectrum and ax!=ax_phase:
            for label in ax.get_xticklabels():
                label.set_visible(False)

    #
    
    
    fig.subplots_adjust(hspace=0)
 
 
#    t_array = np.linspace(g.t[0], g.t[-1], len(g.t))
#    s_array = t_array*c*1.0e-15
#    ncar = g('ncar')
#    dgrid = g('dgrid')
#
#    xlamds = float(g('xlamds'))
#    nslice = len(g.spec)
#    zsep  = float(g('zsep'))
#    zstop  = float(g('zstop'))
#
#    srange = nslice*zsep*xlamds     #range in s
#    dk = 2*np.pi/srange
#    ds = zsep*xlamds
#
#    #w_l_m  =  g('xlamds')
#    w_l_ev = h * c / g('xlamds')

    # max_pw_list = np.zeros(run_end+1)
    # max_argpw_list = np.zeros(run_end+1)
    # max_sp_list = np.zeros(run_end+1)
    # max_argsp_list = np.zeros(run_end+1)

    #

    # plot_lattice(g, font_size)

    # plot_betas(ax_sizepos, S, beta_x, beta_y, font_size)

    # plot_elems(ax_energy, lat, s_point = S[0], legend = legend, y_scale=0.8) # plot elements
    s=g.t*c*1.0e-15*1e6
    
    
        
    
    ax_curr.plot(s, g.I/1e3, 'k--')
    ax_curr.set_ylabel('I[kA]')
    
    ax_time = ax_curr.twiny()

#    X = np.linspace(0,1,1000)
#    Y = np.cos(X*20)
#    
#    ax1.plot(X,Y)
#   #new_tick_locations = np.array([4, 8, 12, 16])
    
    def tick_function(X):
        V = X*3.33
        return ["%.1f" % z for z in V]
    
    ax_time.set_xlim(ax_curr.get_xlim())
#    ax_time.set_xticks(new_tick_locations)
    ax_time.set_xticklabels(tick_function(s))
#    ax_time.set_xticks(new_tick_locations)
    ax_time.set_xlabel(r"t [fs]")
    
    #ax_pow.set_yscale('log')

#    print s.shape 
#    print g.p_int[:,zi].shape

    
#    ax_power.plot = (s, g.p_int[:,zi], 'k-')
    ax_power = ax_curr.twinx()
    ax_power.plot(s,g.p_int[:,zi],'g-',linewidth=1.5)    
    ax_power.set_ylabel('Power [W]')
    ax_power.set_ylim([0, np.amax(g.p_int[:,zi])])
#    ax_power.get_xaxis().get_offset_text().set_x(1.1)

    ax_energy.plot(s, g.el_energy[:,zi]*0.511e-3, 'b-', s, (g.el_energy[:,zi]+g.el_e_spread[:,zi])*0.511e-3, 'r--',s, (g.el_energy[:,zi]-g.el_e_spread[:,zi])*0.511e-3, 'r--',)
    ax_energy.set_ylabel('$E\pm\sigma_E$\n[GeV]')
#    ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.ticklabel_format(useOffset=False, style='plain')   
    phase=unwrap(g.phi_mid[:,zi])
    ax_phase.plot(s, phase-np.roll(phase,1,0), 'k-',linewidth=0.5)
    ax_phase.set_ylabel('$\phi$ [rad]')
    ax_phase.set_ylim([-0.5, 0.5])
    ax_phase.set_xlabel('[$\mum$]')   
    
    n_pad=1
    # print len(g.z),len(g.xrms[0,:]),len(np.mean(g.yrms,axis=0))
    power=np.pad(g.p_mid, [(int(g.nSlices/2)*n_pad, (g.nSlices-(int(g.nSlices/2))))*n_pad, (0, 0)], mode='constant')
    phase=np.pad(g.phi_mid, [(int(g.nSlices/2)*n_pad, (g.nSlices-(int(g.nSlices/2))))*n_pad, (0, 0)], mode='constant')
    # spectrum = np.power(abs(fft(np.sqrt( np.array(g.p_mid)) * np.exp( 1.j* np.array(g.phi_mid) ) , axis=0)),2)#/sqrt(g.nSlices)
    spectrum = np.power(abs(fft(np.sqrt( np.array(power)) * np.exp( 1.j* np.array(phase) ) , axis=0)),2)/sqrt(g.nSlices)/np.power((2*g.leng/g('ncar')),2)/1e10
    e_0=1239.8/g('xlamds')/1e9
    # print e_0

    g.freq_ev1 = h * fftfreq(len(spectrum), d=g('zsep') * g('xlamds') / c)+e_0
    lamdscale=1239.8/g.freq_ev1
    #sys.exit()
    #print 'shape',spectrum.shape

    spectrum_average_z=np.sum(spectrum,axis=0)
    spectrum_average_z[np.where(spectrum_average_z==0)[0]]=1 #avoid division by zero
    spectrum_lamdpos=sum(np.matlib.repmat(lamdscale,spectrum.shape[1],1)*transpose(spectrum),axis=1)/np.sum(spectrum,axis=0)

#    print 'p_mid',g.p_mid.shape
#    print 'lamdpos',spectrum_lamdpos.shape
#    print 'lamdscale',lamdscale.shape    
#    print 'spectrum',spectrum.shape

    ax_spectrum.plot(fftshift(lamdscale), fftshift(spectrum[:,zi]), 'r-')
    ax_spectrum.set_ylabel('P($\lambda$)')
    ax_spectrum.set_xlabel('$\lambda$ [nm]')
    
    #fix!!!
    # spectrum_lamdpos=sum(np.matlib.repmat(lamdscale,spectrum.shape[1],1)*transpose(spectrum),axis=1)/np.sum(spectrum,axis=0)
    # spectrum_width=sqrt(sum((np.power(np.matlib.repmat(lamdscale,spectrum.shape[1],1)-np.matlib.repmat(spectrum_lamdpos,spectrum.shape[1],1),2)*spectrum),1)/sum(spectrum,axis=0));
    # ax_spec_bandw = ax_spectrum.twinx()
    # ax_spec_bandw.plot(g.z, spectrum_lamdpos, 'g-')
    # fix and include!!!



    # print 'spec1', len(spectrum)
    # print 'freq_ev1', len(g.freq_ev1)


    #
    # fig2=plt.plot(fftshift(g.freq_ev1),fftshift(spectrum[:,-1],axes=0))
    # plt.show()
    # fig3=plt.plot(fftshift(lamdscale),fftshift(spectrum[:,-1],axes=0))
    # plt.show()


    #ax_energy.plot(g.z, np.mean(g.energy_GeV,axis=0), 'g-')
    # ax_spread = ax_energy.twinx()
    # ax_spread.plot(g.z, np.average(g.e_spread_GeV, weights=g.I[1:], axis=0), 'b--', g.z, np.amax(g.e_spread_GeV, axis=0), 'b-')
    # ax_spread.set_ylabel('Energy spread ($\Delta$E)')

    number_ticks=4



    # ax_bunching.plot(g.z, np.average(g.bunching, weights=g.I[1:], axis=0), 'k--', g.z, np.amax(g.bunching, axis=0), 'r-')
    # ax_bunching.set_ylabel('Bunching')

    # ax_bunching.set_xlabel('s [m]')
    # x1,x2,y1,y2 = ax_sizepos.axis()
    # ax_sizepos.axis([x1,x2,0,y2])
    # x1,x2,y1,y2 = ax_spread.axis()
    # ax_spread.axis([x1,x2,0,y2])

    # ax_pow
    # ax_spectrum
    # ax_size


    # ax_pow.yaxis.major.locator.set_params(nbins=number_ticks)
    # ax_en.yaxis.major.locator.set_params(nbins=number_ticks)
    # ax_spectrum.yaxis.major.locator.set_params(nbins=number_ticks)
    
    ax_phase.xaxis.major.locator.set_params(nbins=5)
    ax_power.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_spectrum.yaxis.major.locator.set_params(nbins=number_ticks)
    # ax_bunching.yaxis.major.locator.set_params(nbins=number_ticks)

    # yloc = plt.MaxNLocator(max_yticks)
    # # print yloc
    # ax_sizepos.yaxis.set_major_locator(yloc)
    # ax_energy.yaxis.set_major_locator(yloc)
    # ax_und.yaxis.set_major_locator(yloc)
    # ax_bunching.yaxis.set_major_locator(yloc)

    # ax_energy.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1e'))

    plt.xlim(s[0], s[-1])

    fig.subplots_adjust(top=0.95, bottom=0.2, right=0.85, left=0.15)

    #fig.set_size_inches((8,8),forward=True)
    
    pos1 = ax_spectrum.get_position() # get the original position 
    pos2 = [pos1.x0 + 0, pos1.y0 - 0.1,  pos1.width / 1.0, pos1.height / 0.9] 
    ax_spectrum.set_position(pos2)

    ax_spectrum.tick_params(axis='y', which='both', colors='r')
    ax_spectrum.yaxis.label.set_color('r')    
    ax_energy.tick_params(axis='y', which='both', colors='b')
    ax_energy.yaxis.label.set_color('b')    
    ax_power.tick_params(axis='y', which='both', colors='g')
    ax_power.yaxis.label.set_color('g')    
    
    
    
#    ax_spectrum.spines['left'].set_color('red')

    plt.show()
    
#    ax.ticklabel_format(useOffset=False, style='plain')

def plot_lattice(g, font_size):

    plot(srange, g.aw, 'b-')
    xlabel('s [m]')
    ylabel('K_value')

    # turn off the 2nd axes rectangle with frameon kwarg
    ax2 = twinx()
    s2 = sin(2*pi*t)
    plot(t, s2, 'r.')
    ylabel('sin')
    ax2.yaxis.tick_right()
    show()



    ax.set_ylabel(r"$\beta_{x,y}$, m")
    ax.plot(S, beta_x,'r', lw = 2, label=r"$\beta_{x}$")
    ax.plot(S, beta_y,'b', lw = 2, label=r"$\beta_{y}$")
    leg = ax.legend(loc='upper right', shadow=True, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    leg.get_frame().set_alpha(0.5)




    # S = [p.s for p in tws]#map(lambda p:p.s, tws)
    # d_Ftop = []
    # Fmin = []
    # Fmax = []
    # for elem in top_plot:
    #     Ftop = [p.__dict__[elem] for p in tws]# map(lambda p:p.__dict__[elem], tws)
    #     #for f in Ftop:
    #     #    print(f)
    #     #print (max(Ftop))
    #     Fmin.append(min(Ftop))
    #     Fmax.append(max(Ftop))
    #     top_label = r"$"+elem+"$"
    #     ax.plot(S, Ftop, lw = 2, label=top_label)
    #     d_Ftop.append( max(Ftop) - min(Ftop))
    # d_F = max(d_Ftop)
    # if d_F == 0:
    #     d_Dx = 1
    #     ax.set_ylim(( min(Fmin)-d_Dx*0.1, max(Fmax)+d_Dx*0.1))
    # if top_plot[0] == "E":
    #     top_ylabel = r"$"+"/".join(top_plot) +"$"+ ", GeV"
    # else:
    #     top_ylabel = r"$"+"/".join(top_plot) +"$"+ ", m"
    #
    # yticks = ax.get_yticks()
    # yticks = yticks[2::2]
    # ax.set_yticks(yticks)
    # #for i, label in enumerate(ax.get_yticklabels()):
    # #    if i == 0 or i == 1:
    # #        label.set_visible(False)
    # ax.set_ylabel(top_ylabel)
    #
    # #ax.plot(S, Dx,'black', lw = 2, label=lable)
    # leg2 = ax.legend(loc='upper right', shadow=True, fancybox=True,prop=font_manager.FontProperties(size=font_size))
    # leg2.get_frame().set_alpha(0.5)












# import matplotlib.font_manager as font_manager
# font = {
#         'size'   : 20}
# matplotlib.rc('font', **font)


# '''
# try:
    # from PyQt4.QtCore import *
    # from PyQt4.QtGui import *
    # from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    # from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
# except:
    # print 'WARNING: Qt not installed, some graphics may not work properly'
# '''



# def plot_lattice(lat, axis, alpha=1.0, params={'kmax':2.0, 'ang_max':0.5e-2}, s_start = 0.0):
    # axis.grid(True)
    

    # pos = s_start
    # offs = np.array([0,0.0])
     
    # ang_max = params['ang_max']           # max dipole strength in lattice
    # min_dipole_height = 0.1  # dipole of any strength will show at least this strength 

    # min_solenoid_height = 0.1
    # sol_max = 0.1

    # kmax = params['kmax']
    # min_quad_height = 0.1


    # rendered_seq = []
    # rendered_len = 0.0
    # total_len = 0.0
    
    
    # for e in lat.sequence:
        
        # if e.type in ['bend','sbend', 'rbend', 'quadrupole', 'undulator', 'drift', 'monitor','hcor','vcor', 'cavity','edge', 'solenoid']:
            # e.s = total_len
            # rendered_seq.append(e)
            # rendered_len += e.l
        # total_len += e.l
    
    # for e in rendered_seq:
        # dx = e.l 
        
        # if e.type in ['bend', 'sbend', 'rbend','hcor','vcor']:
            # axis.add_patch( mpatches.Rectangle(offs+np.array([pos, 0.0]), dx, 
                                               # np.sign(-e.angle) * min_dipole_height - e.angle / ang_max * (1- min_dipole_height), 
                                               # color='#0099FF', alpha = alpha))

        # if e.type in ['solenoid']:
            # axis.add_patch( mpatches.Rectangle(offs+np.array([pos, 0.0]), dx, 
                                               # np.sign(-e.k) * min_solenoid_height - e.k / sol_max * (1- min_solenoid_height), 
                                               # color='#FF99FF', alpha = alpha))


        # if e.type in ['quadrupole']:
                        
            # if e.k1 >= 0:
                # axis.add_patch( mpatches.Ellipse(offs+np.array([pos, 0.0]), dx, min_quad_height + abs(e.k1/kmax)*2, color='green', alpha=alpha) )
            # else:
                # Path = mpath.Path
                # h = abs(e.k1/kmax) + min_quad_height / 2

                # verts = np.array([
                                  # (dx, h),
                                  # (-dx, h),
                                  # (-dx/4, 0),
                                  # (-dx, -h),
                                  # (dx, -h),
                                  # (dx/4, 0),
                                  # (dx, h)
                                  # ])
                
                # codes = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.LINETO, Path.LINETO, Path.CURVE3, Path.CURVE3]
                
                
                # path = mpath.Path(verts+offs+np.array([pos, 0.0]), codes)
                # patch = mpatches.PathPatch(path, color='green', alpha=alpha)
                # axis.add_patch(patch)

            
        
        # if e.type in ['undulator']:
            # nper = 16
            # dxs = dx / nper / 2.0
            
            # height = 0.7
            # gap = 0.1
            # if e.Kx**2 + e.Ky**2 < 1.e-6:
                # gap = 0.2
            
            # for iseg in np.arange(0,nper,2):                                
                # axis.add_patch( mpatches.Rectangle(offs+np.array([pos + iseg * dxs, gap]), dxs, height, color='red', alpha = alpha))
                # axis.add_patch( mpatches.Rectangle(offs+np.array([pos + iseg * dxs, -height-gap]), dxs, height, color='blue', alpha = alpha) )
                # axis.add_patch( mpatches.Rectangle(offs+np.array([pos + (iseg+1) * dxs, gap]), dxs, height, color='blue', alpha = alpha) )
                # axis.add_patch( mpatches.Rectangle(offs+np.array([pos + (iseg+1) * dxs, -height-gap]), dxs, height, color='red', alpha = alpha  ) )

        # if e.type in ['cavity']:
            # nper = 16
            # dxs = dx / nper / 2.0
            
            # height = 0.7
            # gap = 0.1
            # axis.add_patch( mpatches.Ellipse(offs+np.array([pos + dx/2, 0.0]), 
                            # dx/3, 0.5, 
                            # color='#FF0033', alpha = alpha))

            # axis.add_patch( mpatches.Ellipse(offs+np.array([pos + dx/6, 0.0]), 
                            # dx/3, 0.5, 
                            # color='#FF0033', alpha = alpha))

            # axis.add_patch( mpatches.Ellipse(offs+np.array([pos + 5*dx/6, 0.0]), 
                            # dx/3, 0.5, 
                            # color='#FF0033', alpha = alpha))

     
        # if e.type in ['monitor']:
            # Path = mpath.Path
            # h = 0.005
            # offs_mon = np.array([0.0, 0.03])

            # verts = np.array([
                              # (0, h),
                              # (0, -h),
                              # (h, 0),
                              # (-h, 0)
                              # ])
            
            # codes = [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO]
            
            
            # verts += offs_mon
            
            # path = mpath.Path(verts+offs+np.array([pos, 0.0]), codes)
            # patch = mpatches.PathPatch(path, color='black', lw=2, alpha = alpha*0.5)
            # axis.add_patch(patch)
            # axis.add_patch( mpatches.Circle(offs+offs_mon+np.array([pos, 0.0]), h/2, color='black', alpha = alpha*0.5))

        # pos += dx
    
    
    # lw = 0.005
    # axis.add_patch( mpatches.Rectangle(offs - (0,lw/2), 1.0, lw, color='black', alpha = alpha*0.4) )

 
    # axis.set_ylim(-1,1)
    # axis.set_xticks([])
    # axis.set_yticks([])

# # Tomin Sergey functions


# def elem_cord(lat):
    # quad = np.array([[0,0]])
    # bend = np.array([[0,0]])
    # sext = np.array([[0,0]])
    # multi = np.array([[0,0]])
    # c = []
    # corr = np.array([[0,0]])
    # mons = np.array([[0,0]])
    # cav = np.array([[0,0]])
    # mat = np.array([[0,0]])
    # und = np.array([[0,0]])
    # drft = np.array([[0,0]])
    # L = 0
    # #for elem in lat.sequence:
    # #if elem.type == "drift" and elem.l == 0:
    # #lat.sequence.remove(elem)
    
    # for elem in lat.sequence:
        # dL = 0.
        # if elem.l == 0:
            # dL = 0.03
        # temp = np.array([[L - dL, 0]])
        # temp = np.append(temp, [[L - dL, 1]], axis = 0)
        # temp = np.append(temp, [[L+elem.l + dL, 1]], axis = 0)
        # temp = np.append(temp, [[L+elem.l + dL, 0]], axis = 0)
        # #temp = temp.reshape((4,2))
        
        # if elem.type == "quadrupole":
            # k1 = elem.k1
            # quad = np.append(quad, [[L, 0]], axis=0)
            # quad = np.append(quad, [[L, k1]], axis=0)
            # #quad = np.append(quad, [[L+elem.l/2, k1*(1 + 0.2 * np.sign(elem.k1)) ]], axis=0)
            # quad = np.append(quad, [[L+elem.l, k1]],axis=0)
            # quad = np.append(quad, [[L+elem.l, 0]],axis=0)

        # elif elem.type == "cavity":
            # k1 = 1.
            # cav = np.append(cav, [[L, 0]], axis=0)
            # cav = np.append(cav, [[L, k1]], axis=0)
            # cav = np.append(cav, [[L+elem.l, k1]],axis=0)
            # cav = np.append(cav, [[L+elem.l, 0]],axis=0)
        # elif elem.type == "drift":
            # k1 = 1.
            # drft = np.append(drft, [[L, 0]], axis=0)
            # drft = np.append(drft, [[L, 0]], axis=0)
            # drft = np.append(drft, [[L+elem.l, 0]],axis=0)
            # drft = np.append(drft, [[L+elem.l, 0]],axis=0)
        # elif elem.type == "matrix":
            # #k1 = 1.
            # mat = np.append(mat, [[L, 0]], axis=0)
            # #mat = np.append(mat, [[L, k1]], axis=0)
            # #mat = np.append(mat, [[L+elem.l, k1]],axis=0)
            # mat = np.append(mat, [[L+elem.l, 0]],axis=0)

        # elif elem.type == "undulator":
            # k1 = 1.
            # und = np.append(und, [[L, 0]], axis=0)
            # und = np.append(und, [[L, k1]], axis=0)
            # und = np.append(und, [[L+elem.l, k1]],axis=0)
            # und = np.append(und, [[L+elem.l, 0]],axis=0)

        # elif elem.type in ["sbend", "bend", "rbend"]:
            # if elem.l == 0:
                # h = 0
            # else:
                # h = elem.angle/elem.l
            # temp[:,1] = temp[:,1]*h
            # bend = np.append(bend, temp, axis = 0)
        
        # elif elem.type == "sextupole":

            # temp[:,1] = temp[:,1]*elem.ms+ temp[:,1]*elem.k2
            # sext = np.append(sext, temp, axis = 0)

        # elif elem.type == "multipole":
            # if sum(abs(elem.kn)) != 0:
                # temp[:,1] = temp[:,1]*sum(elem.kn)/sum(abs(elem.kn))
            # else:
                # temp[:,1] = temp[:,1]*0.
            # multi = np.append(multi, temp, axis=0)
        
        # elif elem.type in ["hcor" , "vcor"]:
            # temp[:,1] = temp[:,1]#*abs(elem.angle)
            # corr = np.append(corr, temp, axis=0)
        
        # elif elem.type in ["monitor"]:
            # temp[:,1] = temp[:,1]#*abs(elem.angle)
            # mons = np.append(mons, temp, axis=0)
        # #c.append((L,  elem.l+0.03))
        # L += elem.l
    # if len(quad) != 1:
        # quad[:,1] = quad[:,1]/max(quad[:,1])
    # if len(bend) != 1:
        # if max(bend[:,1]) == 0:
            # bend[:,1] = 0
        # else:
            # bend[:,1] = bend[:,1]/max(bend[:,1])
    # if len(sext) != 1 and max(sext[:,1] != 0):
        # sext[:,1] = sext[:,1]/max(np.abs(sext[:,1]))
    # if len(corr) != 1 and max(corr[:,1] != 0):
        # corr[:,1] = corr[:,1]/max(corr[:,1])
    # #if len(corr) != 1 and max(mons[:,1] != 0):
    # #mons[:,1] = mons[:,1]/max(mons[:,1])
    # return quad, bend, sext, corr, mons, cav, mat, und, multi, drft


# def plot_elems(ax, lat, s_point = 0, nturns = 1, y_lim = None,y_scale = 1, legend = True):
    # quad, bend, sext, corr, mons, cav, mat, und, multi, drft = elem_cord(lat)
    # #print len(quad), len(bend), len(sext), len(corr ),len( mons), len( cav)
    # #print cav
    # alpha = 1
    # ax.set_ylim((-1,1.5))
    # if y_lim != None:
        # ax.set_ylim(y_lim)
    # n = 0
    # for i in range(nturns):
        # n = 0
        # if len(quad)>1:
            # ax.fill(quad[:,0]+i*lat.totalLen + s_point, quad[:,1]*y_scale*0.8, "r", alpha = alpha, label = "quad")
            # n += 1
        # if len(bend)>1:
            # ax.fill(bend[:,0]+i*lat.totalLen + s_point, bend[:,1]*y_scale*0.7, "lightskyblue", alpha = alpha, label = "bend")
            # n += 1
        # if len(sext)>1:
            # ax.fill(sext[:,0]+i*lat.totalLen + s_point, sext[:,1]*y_scale*0.8, "green",  edgecolor = "green", alpha = alpha, label = "sext")
            # n += 1
        # if len(multi)>1:
            # ax.fill(multi[:,0]+i*lat.totalLen + s_point, multi[:,1]*y_scale*0.8, "green",  edgecolor = "green", alpha = alpha, label = "multi")
            # n += 1
        # if len(corr)>1:
            # ax.fill(corr[:,0]+i*lat.totalLen + s_point, corr[:,1]*y_scale*0.7, "b", edgecolor = "b", alpha = alpha, label = "corr")
            # n += 1
        # if len(mons)>1:
            # ax.fill(mons[:,0]+i*lat.totalLen + s_point, mons[:,1]*y_scale*0.7, "orange", edgecolor = "orange", alpha = alpha, label = "bpm")
            # n += 1
        # if len(cav)>1:
            # ax.fill(cav[:,0]+i*lat.totalLen + s_point, cav[:,1]*y_scale*0.7, "orange", edgecolor = "lightgreen", alpha = alpha, label = "cav")
            # #print cav[:,0]+i*lat.totalLen, cav[:,1]*y_scale*0.7, i*lat.totalLen
            # n += 1
        # if len(mat)>1:
            # ax.fill(mat[:,0]+i*lat.totalLen + s_point, mat[:,1]*y_scale*0.7, "pink", edgecolor = "lightgreen", alpha = alpha, label = "mat")
            # n += 1
        # if len(und)>1:
            # ax.fill(und[:,0]+i*lat.totalLen + s_point, und[:,1]*y_scale*0.7, "pink", edgecolor = "lightgreen", alpha = alpha, label = "und")
            # n += 1
        # if len(drft)>1:
            # ax.fill(drft[:,0]+i*lat.totalLen + s_point, drft[:,1]*y_scale*0.7, "k")
            # n += 1
    # #ax.broken_barh(s , (y0, yw*1.3), facecolors='green', edgecolor = "green", alpha = alpha, label = "Sext")
    # #ax.broken_barh(c , (y0, yw), facecolors='b',edgecolor = "b", alpha = alpha)
    # if legend:
        # ax.legend(loc='upper center', ncol=n, shadow=False, prop=font_manager.FontProperties(size=15))

# def plot_disp(ax,tws, top_plot, font_size):
    # S = [p.s for p in tws]#map(lambda p:p.s, tws)
    # d_Ftop = []
    # Fmin = []
    # Fmax = []
    # for elem in top_plot:
        # Ftop = [p.__dict__[elem] for p in tws]# map(lambda p:p.__dict__[elem], tws)
        # #for f in Ftop:
        # #    print(f)
        # #print (max(Ftop))
        # Fmin.append(min(Ftop))
        # Fmax.append(max(Ftop))
        # top_label = r"$"+elem+"$"
        # ax.plot(S, Ftop, lw = 2, label=top_label)
        # d_Ftop.append( max(Ftop) - min(Ftop))
    # d_F = max(d_Ftop)
    # if d_F == 0:
        # d_Dx = 1
        # ax.set_ylim(( min(Fmin)-d_Dx*0.1, max(Fmax)+d_Dx*0.1))
    # if top_plot[0] == "E":
        # top_ylabel = r"$"+"/".join(top_plot) +"$"+ ", GeV"
    # else:
        # top_ylabel = r"$"+"/".join(top_plot) +"$"+ ", m"

    # yticks = ax.get_yticks()
    # yticks = yticks[2::2]
    # ax.set_yticks(yticks)
    # #for i, label in enumerate(ax.get_yticklabels()):
    # #    if i == 0 or i == 1:
    # #        label.set_visible(False)
    # ax.set_ylabel(top_ylabel)
    
    # #ax.plot(S, Dx,'black', lw = 2, label=lable)
    # leg2 = ax.legend(loc='upper right', shadow=True, fancybox=True,prop=font_manager.FontProperties(size=font_size))
    # leg2.get_frame().set_alpha(0.5)




# def plot_betas(ax, S, beta_x, beta_y, font_size):
    # ax.set_ylabel(r"$\beta_{x,y}$, m")
    # ax.plot(S, beta_x,'r', lw = 2, label=r"$\beta_{x}$")
    # ax.plot(S, beta_y,'b', lw = 2, label=r"$\beta_{y}$")
    # leg = ax.legend(loc='upper right', shadow=True, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    # leg.get_frame().set_alpha(0.5)


# def plot_xy(ax, S, X, Y, font_size):
    # ax.set_ylabel(r"$X, Y$, m")
    # ax.plot(S, X,'r', lw = 2, label=r"$X$")
    # ax.plot(S, Y,'b', lw = 2, label=r"$Y$")
    # leg = ax.legend(loc='upper right', shadow=True, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    # leg.get_frame().set_alpha(0.5)




# def body_trajectory(fig, ax_xy, ax_el, lat, list_particles):
    # X = map(lambda p:p.x, list_particles)
    # Y = map(lambda p:p.y, list_particles)
    # S = map(lambda p:p.s, list_particles)
    
    # font_size = 16
    
    # for ax in ax_xy, ax_el:
        # if ax!=ax_el:
            # for label in ax.get_xticklabels():
                # label.set_visible(False)
    
    # ax_xy.grid(True)
    # ax_el.set_yticks([])
    # ax_el.grid(True)
    # plt.xlim(S[0], S[-1])
    
    # fig.subplots_adjust(hspace=0)
    
    # plot_xy(ax_xy, S, X, Y, font_size)

    # plot_elems(ax_el, lat, nturns = int(S[-1]/lat.totalLen), legend = False) # plot elements

# """
# def plot_current(p_array, charge, num_bins = 200):
    # z = p_array.particles[4::6]
    # hist, bin_edges = np.histogram(z, bins=num_bins)
    # delta_Z = max(z) - min(z)
    # delta_z = delta_Z/num_bins
    # t_bins = delta_z/speed_of_light
    # print "Imax = ", max(hist)*charge/t_bins
    # hist = np.append(hist, hist[-1])
    # plt.plot(bin_edges, hist*charge/t_bins)
    # plt.grid(True)
    # plt.title("current")
    # plt.xlabel("s, m")
    # plt.ylabel("I, A")
    # plt.show()
# """

# def plot_trajectory(lat, list_particles):
    # fig = plt.figure()
    # plt.rc('axes', grid=True)
    # plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    # left, width = 0.1, 0.85
    # rect2 = [left, 0.2, width, 0.7]
    # rect3 = [left, 0.05, width, 0.15]
    
    # ax_xy = fig.add_axes(rect2)  #left, bottom, width, height
    # ax_el = fig.add_axes(rect3, sharex=ax_xy)
    
    # body_trajectory(fig, ax_xy, ax_el, lat, list_particles)
    # plt.show()

# def plot_trajectory_test(fig, lat, p1, p2, p3,p4, alpha = 1):
    # #fig = plt.figure()
    # plt.rc('axes', grid=True)
    # plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    # left, width = 0.1, 0.85
    # rect2 = [left, 0.2, width, 0.7]
    # rect3 = [left, 0.05, width, 0.15]
    
    # ax_xy = fig.add_axes(rect2)  #left, bottom, width, height
    # ax_el = fig.add_axes(rect3, sharex=ax_xy)
    
    # X = array(map(lambda p:p.x, p1))
    # S = map(lambda p:p.s, p1)
    # X2 = array(map(lambda p:p.x, p2))
    # S2 = map(lambda p:p.s, p2)
    
    # font_size = 16
    
    # for ax in ax_xy, ax_el:
        # if ax!=ax_el:
            # for label in ax.get_xticklabels():
                # label.set_visible(False)
    
    # ax_xy.grid(True)
    # ax_el.set_yticks([])
    # ax_el.grid(True)
    # plt.xlim(S[0], S[-1])
    
    # fig.subplots_adjust(hspace=0)
    # Si = map(lambda p:p.s, p3)
    # Xi1 = array(map(lambda p:p.x, p3))*1000
    # Xi2 = array(map(lambda p:p.x, p4))*1000
    # ax_xy.plot(Si, Xi1,'b', lw = 2, label=r"inj $\pm \sigma_x$")
    # ax_xy.plot(Si, Xi2,'b', lw = 2)
    # ax_xy.fill_between(Si, Xi1,Xi2, alpha=alpha, facecolor='blue', label =r"inj $\pm  \sigma_x$")
    # #plot_xy(ax_xy, S, X, Y, font_size)
    # ax_xy.set_ylabel(r"$X$, mm")
    
    # ax_xy.plot(S, X*1000,'r', lw = 2, label=r"store $\pm \sigma_x$")
    # ax_xy.plot(S2, X2*1000,'r', lw = 2)
    # ax_xy.fill_between(S, X*1000,X2*1000, alpha=alpha, facecolor='red', label =r"store $\pm  \sigma_x$")
    
    
    
    
    # ax_xy.broken_barh([(10, 0.344)] , (-20, -2.4), facecolors='black')
    # ax_xy.broken_barh([(10, 0.344)] , (-22.4, -10), facecolors='yellow')
    # #ax_xy.plot([10., 10, 10.344, 10.344, 10.], array([0.02,0.0224, 0.0224, 0.02, 0.02])*1000, 'k')
    # #ax_xy.plot([10., 10, 10.344, 10.344, 10.], array([0.0224,0.0324, 0.0324, 0.0224, 0.0224])*1000, 'k')
    # leg = ax_xy.legend(loc='upper right', shadow=True, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    # leg.get_frame().set_alpha(0.5)
    # plot_elems(ax_el, lat, nturns = int(S[-1]/lat.totalLen), legend = True) # plot elements
# #plt.show()

# def plot_traj_pulse(lat, list_particles, list_particles2, U1, U2):
    # fig = plt.figure()
    # plt.rc('axes', grid=True)
    # plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    # left, width = 0.1, 0.85
    # rect2 = [left, 0.2, width, 0.7]
    # rect3 = [left, 0.05, width, 0.15]
    
    # ax_xy = fig.add_axes(rect2)  #left, bottom, width, height
    # ax_el = fig.add_axes(rect3, sharex=ax_xy)
    # S = map(lambda p:p.s, list_particles)
    # ax_xy.plot(S, U1,'g', lw = 1, label=r"$U1$")
    # ax_xy.plot(S, U2,'b', lw = 1, label=r"$U2$")
    # #print map(lambda p:p.s, list_particles2)
    # ax_xy.plot(map(lambda p:p.s, list_particles2), map(lambda p:p.x, list_particles2),'b', lw = 2, label=r"$p2$")
    # body_trajectory(fig, ax_xy, ax_energy, lat, list_particles)
    # #body_trajectory(fig, ax_xy, ax_energy, lat, list_particles2)
    # plt.show()



# def plot_elem_disp(lat, tws):
    # Dx = map(lambda p:p.Dx, tws)
    # S = map(lambda p:p.s, tws)
    # fig, ax = plt.subplots(1,sharex=True)
    # ax.grid(False)
    # ax.set_yticks([])
    # plt.xlim(tws[0].s, tws[-1].s)
    # lim = max(Dx)*1.1
    # plot_disp(ax, S, Dx, font_size = 16)
    # plot_elems(ax, lat, y_lim = None, y_scale = lim*0.6, legend = False) # plot elements
    # ax.set_ylim((0,lim))
    # plt.show()


# def resonans(Qx, Qy, order = 5):
    # ORD = order
    # qx1, qy1 = 0,0
    # qx2, qy2 = 2,2
    # X = []
    # Y = []
    # Order = []
    # params = []
    # #Qx = 0.22534
    # #Qy = 0.301
    # for order in range(1, ORD + 1):
        # n = np.arange(-order, order+1)
        # m = np.array([order - abs(n), - order + abs(n)]).flatten()
        # n = np.tile(n, 2)
        # ziped = []
        # for n,m in zip(n,m):
            # if (n,m) in ziped:
                # continue
            # ziped.append((n,m))
        # #ziped =  zip(n,m)
        # #print ziped
        # for p in range(-50, 50):
            # for n, m in ziped:
                # #print p, n,m
                # if m != 0:
                    # x = [qx1, qx2]
                    # y = [(p - n*qx1)/float(m), (p - n*qx2)/float(m)]
                # else:
                    # x = [p/float(n), p/float(n)]
                    # y = [qy1, qy2]
                # params.append([n,m,p])
                # X.append(x)
                # Y.append(y)
                # Order.append(order)
    # return X,Y,Order,params

# def resonans_diag(Qx, Qy, order):
    # X,Y,Order,params = resonans(Qx, Qy, order)
    # indsort = np.argsort(Order)
    # #print Order
    # #print len(indsort), len(X)
    # X = np.array(X)
    # #print X[indsort]
    # #fig = plt.figure()
    # #ax1 = fig.add_subplot(111)
    # """
        # def onpick3(event):
        # ind = event.ind
        # print 'onpick3 scatter:', ind #, np.take(x, ind), np.take(y, ind)
        
        # ax1.plot(X[:], Y[:], picker=True)
        # plt.xlim(0,1)
        # plt.ylim(0,1)
        # """
    # for i, order in enumerate(Order):
        # if order == 1:
            # color = "k"
            # lw = 3
        # elif order == 2:
            # color ="r"
            # lw = 2
        # elif order == 3:
            # color = "b"
            # lw = 2
        # elif order == 4:
            # color = "g"
            # lw = 0.6
        # else:# order == 5:
            # color = "c"
            # lw = 0.3
        # #print array(X[i])+Qx
        # #print array([i])+Qy
        # plt.plot(array(X[i])+Qx, array(Y[i])+Qy, color, lw = lw, picker=True)
        # plt.xlim(Qx,Qx+1)
        # plt.ylim(Qy,Qy+1)
        # #plt.xticks(x, labels, rotation='vertical')


# def show_da(out_da, x_array, y_array, title=""):
    # from matplotlib import pyplot as plt
    # from numpy import linspace, max, min
    # #print "time execution = ", time() - start , " s"
    # nx = len(x_array)
    # ny = len(y_array)
    # #print(nx, ny, len(out_da))
    # out_da = out_da.reshape(ny,nx)
    # xmin, xmax, ymin, ymax = min(x_array), max(x_array), min(y_array), max(y_array)
    # #plt.subplot(111, axisbg='darkslategray')
    # extent = xmin, xmax, ymin, ymax
    # #print extent
    # #plt.savetxt("da.txt", da)
    # plt.figure(figsize=(10, 7))
    # fig1 = plt.contour(out_da, linewidths=2,extent = extent)#, colors = 'r')
    # #fig1 = plt.contourf(out_da, 20,cmap=plt.cm.rainbow,extent = extent)#, colors = 'r')
    # #plt.axis_bgcolor("#bdb76b")
    # plt.grid(True)
    # plt.title(title)
    # plt.xlabel("X, m")
    # plt.ylabel("Y, m")
    # cb = plt.colorbar()
    # cb.set_label('Nturns')
    # #cb.ax.set_yticklabels(map(str, linspace(min(out_da), max(out_da), 5) ))

    # #plt.savefig('da_error_'+str(int(np.random.rand()*100))+'.png')
    # plt.show()

# def show_mu(contour_da, mux, muy, x_array, y_array, zones = None ):
    # from matplotlib import pyplot as plt

    # nx = len(x_array)
    # ny = len(y_array)
    # t= linspace(0,3.14, num = 100)
    # contour_da = contour_da.reshape(ny,nx)
    # mux = mux.reshape(ny,nx)
    # muy = muy.reshape(ny,nx)
    # xmin, xmax, ymin, ymax = min(x_array), max(x_array), min(y_array), max(y_array)
    # plt.figure(1,figsize=(10, 7)) #axisbg='darkslategray'
    # extent = xmin, xmax, ymin, ymax

    # my_cmap = plt.cm.Paired
    # #my_cmap.set_under('w')
    # #norm = mlb.colors.Normalize(vmin=-0.005, vmax=max(mux))
    # fig1 = plt.contour(contour_da, 1,extent = extent, linewidths=2,colors='k')#, colors = 'r')
    # fig1 = plt.contourf(mux,40, cmap=my_cmap, extent = extent)#, colors = 'r')
    # cb = plt.colorbar(cmap=my_cmap)
    # fig1 = plt.contourf(mux,10, levels=[-1,-.0001], colors='w',extent = extent)
    # if zones != None:
        # x_zone = zones[0]
        # y_zone = zones[1]
        # plt.plot(x_zone*cos(t), y_zone*sin(t), "g", lw = 2)
        # plt.plot(2*x_zone*cos(t), 2*y_zone*sin(t), "b", lw = 2)
        # plt.plot(3*x_zone*cos(t), 3*y_zone*sin(t), "r", lw = 2)
        # plt.plot(4*x_zone*cos(t), 4*y_zone*sin(t), "y", lw = 2)
    # plt.grid(True)
    # #plt.figure(figsize=(10, 7))
    # plt.xlabel("X, m")
    # plt.ylabel("Y, m")
    # cb.set_label('Qx')
    # plt.figure(2,figsize=(10, 7))

    # fig1 = plt.contour(contour_da, 1,extent = extent, linewidths=2,colors='k')#, colors = 'r')
    # fig1 = plt.contourf(muy,40, cmap=my_cmap, extent = extent)#, colors = 'r')
    # if zones != None:
        # x_zone = zones[0]
        # y_zone = zones[1]
        # plt.plot(x_zone*cos(t), y_zone*sin(t), "g", lw = 2)
        # plt.plot(2*x_zone*cos(t), 2*y_zone*sin(t), "b", lw = 2)
        # plt.plot(3*x_zone*cos(t), 3*y_zone*sin(t), "r", lw = 2)
        # plt.plot(4*x_zone*cos(t), 4*y_zone*sin(t), "y", lw = 2)
    # #x = np.linspace(-, 0.01, 0.0001)
    # #plt.plot()
    # cb = plt.colorbar(cmap=my_cmap)
    # fig1 = plt.contourf(muy,10, levels=[-1,-.0001], colors='w',extent = extent)
    # plt.xlabel("X, m")
    # plt.ylabel("Y, m")
    # plt.grid(True)
    # cb.set_label('Qy')
    # plt.show()
