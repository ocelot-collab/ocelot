'''
user interface for viewing genesis simulation results
'''

import sys, os, csv

 
import matplotlib
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from numpy import matlib

from pylab import * #tmp

# font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 20}
params = {'backend': 'ps', 'axes.labelsize': 15, 'font.size': 16, 'legend.fontsize': 24, 'xtick.labelsize': 19,  'ytick.labelsize': 19, 'text.usetex': True}
rcParams.update(params)
rc('text', usetex=True) # required to have greek fonts on redhat

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)

h = 4.135667516e-15
c = 299792458.0

max_yticks = 7

def gen_outplot_e(g, figsize=(8,10), legend = True, fig_name = None, save=False):
    import matplotlib.ticker as ticker

    print('    plotting e-beam evolution')

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

    ax_und=fig.add_subplot(4, 1, 1)
    ax_und.clear()
    ax_size_tpos=fig.add_subplot(4, 1, 2,sharex=ax_und)
    ax_size_tpos.clear()
    ax_energy=fig.add_subplot(4, 1, 3,sharex=ax_und)
    ax_energy.clear()
    ax_bunching=fig.add_subplot(4, 1, 4,sharex=ax_und)
    ax_bunching.clear()

    for ax in ax_size_tpos, ax_energy, ax_und, ax_bunching:
        if ax!=ax_bunching:
            for label in ax.get_xticklabels():
                label.set_visible(False)

    # for tick in ax.yaxis.get_major_ticks():
    #     tick.label.set_fontsize(14)
    #     # specify integer or one of preset strings, e.g.
    #     #tick.label.set_fontsize('x-small')
    #     tick.label.set_rotation('vertical')
    fig.subplots_adjust(hspace=0)

    ax_und.plot(g.z, g.aw, 'b-',linewidth=1.5)
    ax_und.set_ylabel('K (rms)')

    ax_quad = ax_und.twinx()
    ax_quad.plot(g.z, g.qfld, 'r-',linewidth=1.5)
    ax_quad.set_ylabel('Quad')
    ax_quad.grid(False)

    #sys.exit()
    ax_size_tpos.plot(g.z, np.mean(g.xrms,axis=0)*1e6, 'g-',g.z, np.mean(g.yrms,axis=0)*1e6, 'b-')
    ax_size_tpos.set_ylabel('$\sigma_{x,y}$ [$\mu$m]')


    # ax_energy.plot(g.z, np.average(g.el_energy*0.511e-3, weights=g.I, axis=0), 'b-',linewidth=1.5) #with current as weight 
    ax_energy.plot(g.z, np.average(g.el_energy*0.511e-3, axis=0), 'b-',linewidth=1.5)
    ax_energy.set_ylabel('E [GeV]')
    ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_spread = ax_energy.twinx()
    ax_spread.plot(g.z, np.average(g.el_e_spread*0.511e-3*1000, weights=g.I, axis=0), 'm--', g.z, np.amax(g.el_e_spread*0.511e-3*1000, axis=0), 'r--',linewidth=1.5)
    ax_spread.set_ylabel('$\sigma_E$ [MeV]')
    ax_spread.grid(False)
    
    ax_bunching.plot(g.z, np.average(g.bunching, weights=g.I, axis=0), 'k-', g.z, np.amax(g.bunching, axis=0), 'grey',linewidth=1.5)
    ax_bunching.set_ylabel('Bunching')

    ax_bunching.set_xlabel('z [m]')
    
    ax_size_tpos.set_ylim(ymin=0)
    ax_spread.set_ylim(ymin=0)
    ax_bunching.set_ylim(ymin=0)
    
    number_ticks=6
    
    ax_und.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_quad.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_spread.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_bunching.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_size_tpos.yaxis.major.locator.set_params(nbins=number_ticks)
    # yloc = plt.MaxNLocator(max_yticks)
    # ax_size_tpos.yaxis.set_major_locator(yloc)
    # ax_energy.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1e'))
    
    plt.xlim(g.z[0], g.z[-1])
    
    fig.subplots_adjust(top=0.95, bottom=0.1, right=0.85, left=0.15)
    
    aw_tmp=np.array(g.aw)[np.array(g.aw)!=0]
    if np.amax(aw_tmp)!=np.amin(aw_tmp):
        ax_und.set_ylim([np.amin(aw_tmp),np.amax(aw_tmp)])
    ax_und.tick_params(axis='y', which='both', colors='b')
    ax_und.yaxis.label.set_color('b')
    ax_quad.tick_params(axis='y', which='both', colors='r')
    ax_quad.yaxis.label.set_color('r') 
    ax_energy.tick_params(axis='y', which='both', colors='b')
    ax_energy.yaxis.label.set_color('b')    
    ax_spread.tick_params(axis='y', which='both', colors='r')
    ax_spread.yaxis.label.set_color('r') 

    if save!=False:
        if save==True:
            save='png'
        fig.savefig(g.path+'_elec.'+str(save),format=save)
        
    return fig



def gen_outplot_ph(g, figsize=(8, 10), legend = True, fig_name = None, save=False):
    import matplotlib.ticker as ticker
    
    print('    plotting radiation evolution')
    
    font_size = 1
    if fig_name is None:
        if g.filename is '':
            fig = plt.figure('Radaition')
        else:
            fig = plt.figure('Radiation '+g.filename)
    else:
        fig = plt.figure(fig_name)

    fig.set_size_inches(figsize,forward=True)

    plt.rc('axes', grid=True)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    plt.clf()

    
    if g('itdp')==True:
        ax_pow=fig.add_subplot(3, 1, 1)
        ax_pow.clear()
        ax_spectrum=fig.add_subplot(3, 1, 2,sharex=ax_pow)
        ax_spectrum.clear()
        ax_size_t=fig.add_subplot(3, 1, 3,sharex=ax_pow)
        ax_size_t.clear()
        for ax in ax_pow, ax_spectrum, ax_size_t:
            if ax!=ax_size_t:
                for label in ax.get_xticklabels():
                    label.set_visible(False)
    else:
        ax_pow=fig.add_subplot(2, 1, 1)
        ax_pow.clear()
        ax_size_t=fig.add_subplot(2, 1, 2,sharex=ax_pow)
        ax_size_t.clear()
        for ax in ax_pow, ax_size_t:
            if ax!=ax_size_t:
                for label in ax.get_xticklabels():
                    label.set_visible(False)



    # for tick in ax.yaxis.get_major_ticks():
    #     tick.label.set_fontsize(14)
    #     # specify integer or one of preset strings, e.g.
    #     #tick.label.set_fontsize('x-small')
    #     tick.label.set_rotation('vertical')

    #
    fig.subplots_adjust(hspace=0)
        
    ax_pow.plot(g.z, np.amax(g.p_int, axis=0), 'g-',linewidth=1.5)
    ax_pow.set_ylabel('P [W]')
    ax_pow.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_pow.get_yaxis().get_major_formatter().set_scientific(True)
#    if np.amin(g.p_int)>0:
    ax_pow.set_yscale('log')



    ax_en = ax_pow.twinx()
    ax_en.plot(g.z, np.mean(g.p_int,axis=0)*g('xlamds')*g('zsep')*g.nSlices/c, 'k--',linewidth=1.5)
    ax_en.set_ylabel('E [J]')
    ax_en.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_en.get_yaxis().get_major_formatter().set_scientific(True)
#    if np.amin(g.p_int)>0:
    ax_en.set_yscale('log')


    if g('itdp')==True:
        n_pad=1
        # print len(g.z),len(g.xrms[0,:]),len(np.mean(g.yrms,axis=0))
        power=np.pad(g.p_mid, [(int(g.nSlices/2)*n_pad, (g.nSlices-(int(g.nSlices/2))))*n_pad, (0, 0)], mode='constant')
        phase=np.pad(g.phi_mid, [(int(g.nSlices/2)*n_pad, (g.nSlices-(int(g.nSlices/2))))*n_pad, (0, 0)], mode='constant')
        spectrum = abs(fft(np.sqrt( np.array(power)) * np.exp( 1.j* np.array(phase) ) , axis=0))**2/sqrt(g.nSlices)/(2*g.leng/g('ncar'))**2/1e10
        e_0=1239.8/g('xlamds')/1e9
        # print e_0
    
        g.freq_ev1 = h * fftfreq(len(spectrum), d=g('zsep') * g('xlamds') / c)+e_0
        lamdscale=1239.8/g.freq_ev1
        lamdscale_array=np.swapaxes(np.tile(lamdscale,(g.nZ,1)),0,1)    
        
    #    print spectrum.shape
        spectrum_norm=np.sum(spectrum,axis=0)#avoiding division by zero
        spectrum_norm[spectrum_norm==0]=1
    #    print spectrum_norm.shape
        spectrum_lamdpos=np.sum(spectrum*lamdscale_array/spectrum_norm,axis=0)
    #    print "spectrum lamdpos", spectrum_lamdpos
        spectrum_lamdwidth=sqrt(np.sum(spectrum*(lamdscale_array-spectrum_lamdpos)**2/spectrum_norm,axis=0))    
        
        spectrum_lamdwidth1=np.empty(g.nZ)
        for zz in range(g.nZ):
            if np.sum(spectrum[:,zz])!=0:
                peak=fwhm3(spectrum[:,zz])
                #spectrum_lamdwidth1[zz]=abs(lamdscale[peak[0]]-lamdscale[peak[0]+1])*peak[1] #the FWHM of spectral line (error when paekpos is at the edge of lamdscale)
                spectrum_lamdwidth1[zz]=abs(lamdscale[0]-lamdscale[1])*peak[1] #the FWHM of spectral line (error when paekpos is at the edge of lamdscale)
            else:
                spectrum_lamdwidth1[zz]=0
    
    
    
        ax_spectrum.plot(g.z, np.amax(spectrum,axis=0), 'r-',linewidth=1.5)
        ax_spectrum.set_ylabel('P$(\lambda)_{max}$ [a.u.]')
    #    if np.amin(np.amax(spectrum,axis=0))>0:
        ax_spectrum.set_yscale('log')
        
        #fix!!!
        ax_spec_bandw = ax_spectrum.twinx()
        ax_spec_bandw.plot(g.z, spectrum_lamdwidth*2, 'm--')
        ax_spec_bandw.set_ylabel('$2\sigma\lambda$ [nm]')
        # fix and include!!!
    
    
        s=g.t*c*1.0e-15*1e6
        s_array=np.swapaxes(np.tile(s,(g.nZ,1)),0,1)    
        p_int_norm=np.sum(g.p_int,axis=0)#avoiding division by zero
        p_int_norm[p_int_norm==0]=1
        rad_longit_pos=np.sum(g.p_int*s_array/p_int_norm,axis=0)
        rad_longit_size=sqrt(np.sum(g.p_int*(s_array-rad_longit_pos)**2/p_int_norm,axis=0)) #this is standard deviation (sigma)
    
        #g.p_int=np.amax(g.p_int)/1e6+g.p_int # nasty fix from division by zero
        weight=g.p_int+np.amin(g.p_int[g.p_int!=0])/1e6
        
        ax_size_l = ax_size_t.twinx() #longitudinal size
        ax_size_l.plot(g.z, rad_longit_size*2, color='indigo', linestyle='dashed',linewidth=1.5)
        ax_size_l.set_ylabel('longitudinal [$\mu$m]')

        ax_size_t.plot(g.z, np.average(g.r_size*2*1e6, weights=weight, axis=0), 'b-',linewidth=1.5)
        ax_size_t.plot([np.amin(g.z), np.amax(g.z)],[g.leng*1e6, g.leng*1e6], 'b-',linewidth=1.0)
        ax_size_t.set_ylabel('transverse [$\mu$m]')
    else:
        ax_size_t.plot(g.z, g.r_size.T*2*1e6, 'b-',linewidth=1.5)
        ax_size_t.plot([np.amin(g.z), np.amax(g.z)],[g.leng*1e6, g.leng*1e6], 'b-',linewidth=1.0)
        ax_size_t.set_ylabel('transverse [$\mu$m]')



    plt.xlim(g.z[0], g.z[-1])

    fig.subplots_adjust(top=0.95, bottom=0.1, right=0.85, left=0.15)


    ax_pow.tick_params(axis='y', which='both', colors='g')
    ax_pow.yaxis.label.set_color('g')  
    ax_en.tick_params(axis='y', which='both', colors='k')
    ax_en.yaxis.label.set_color('k') 
    ax_en.grid(False)
    ax_size_t.tick_params(axis='y', which='both', colors='b')
    ax_size_t.yaxis.label.set_color('b') 
    ax_size_t.set_xlabel('z [m]')
    ax_size_t.set_ylim(ymin=0)
    ax_pow.yaxis.get_offset_text().set_color(ax_pow.yaxis.label.get_color())
    ax_en.yaxis.get_offset_text().set_color(ax_en.yaxis.label.get_color())

    if g('itdp')==True:
        ax_spectrum.tick_params(axis='y', which='both', colors='r')
        ax_spectrum.yaxis.label.set_color('r')  
        ax_spec_bandw.tick_params(axis='y', which='both', colors='m')
        ax_spec_bandw.yaxis.label.set_color('m')
        ax_spec_bandw.grid(False)
        ax_size_l.tick_params(axis='y', which='both', colors='indigo')
        ax_size_l.yaxis.label.set_color('indigo') 
        ax_size_l.grid(False)
        ax_size_l.set_ylim(ymin=0)    
        ax_spec_bandw.set_ylim(ymin=0)


    # #attempt to fix overlapping label values
# #    for a in [ax_size_l,ax_size_t,ax_spec_bandw,ax_spectrum]:
# #        xticks = a.yaxis.get_major_ticks()
# #        xticks[-1].label.set_visible(False)  
    
# #    labels = ax_size_t.get_yticklabels()
# ##    print dir(labels), labels
# #    labels[0] = ""
# #    ax_size_t.set_yticklabels(labels)    

    if save!=False:
        if save==True:
            save='png'
        fig.savefig(g.path+'_rad.'+str(save),format=save)
    return fig



def gen_outplot_z(g, figsize=(8, 10), legend = True, fig_name = None, z=inf, save=False):
#    max_yticks = 7
    if g('itdp')==False:
        print('    plotting bunch profile at '+str(z)+' [m]')
        print('!     not applicable for steady-state')
        return
    
    import matplotlib.ticker as ticker
    
    if z==inf:
#        print 'Showing profile parameters at the end of undulator'
        z=np.amax(g.z)

    elif z>np.amax(g.z):
#        print 'Z parameter too large, setting to the undulator end'
        z=np.amax(g.z)
    
    elif z<np.amin(g.z):
#        print 'Z parameter too small, setting to the undulator entrance'    
        z=np.amin(g.z)
        
    zi=np.where(g.z>=z)[0][0]
    z=g.z[zi];
    
    print('    plotting bunch profile at '+str(z)+' [m]')



    font_size = 1
    if fig_name is None:
        if g.filename is '':
            fig = plt.figure('Bunch profile at '+str(z)+'m')
        else:
            fig = plt.figure('Bunch profile at '+str(z)+'m '+g.filename)
    else:
        fig = plt.figure(fig_name)
    fig.set_size_inches(figsize,forward=True)
    plt.rc('axes', grid=True)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    # left, width = 0.1, 0.85
    
    plt.clf()
    
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
 
    s=g.t*c*1.0e-15*1e6
    
    
        
    
    ax_curr.plot(s, g.I/1e3, 'k--')
    ax_curr.set_ylabel('I [kA]')
    ax_curr.set_ylim(ymin=0)


    ax_power = ax_curr.twinx()
    ax_power.grid(False)
    ax_power.plot(s,g.p_int[:,zi],'g-',linewidth=1.5)    
    ax_power.set_ylabel('Power [W]')
    ax_power.set_ylim(ymin=0)
    # if np.amax(g.p_int[:,zi])!=np.amin(g.p_int[:,zi]):
        # ax_power.set_ylim([0, np.amax(g.p_int[:,zi])])
    ax_power.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_power.get_yaxis().get_major_formatter().set_scientific(True)
    ax_power.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))#[:,75,75]
    
#    ax_power.get_xaxis().get_offset_text().set_x(1.1)

    ax_energy.plot(s, g.el_energy[:,zi]*0.511e-3, 'b-', s, (g.el_energy[:,zi]+g.el_e_spread[:,zi])*0.511e-3, 'r--',s, (g.el_energy[:,zi]-g.el_e_spread[:,zi])*0.511e-3, 'r--')
    ax_energy.set_ylabel('$E\pm\sigma_E$ [GeV]')
#    ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.ticklabel_format(useOffset=False, style='plain')   

    ax_bunching = ax_energy.twinx()
    ax_bunching.plot(s,g.bunching[:,zi],'grey',linewidth=0.5)
    ax_bunching.set_ylabel('Bunching')
    ax_bunching.set_ylim(ymin=0)
    ax_bunching.grid(False)
    
    
    n_pad=1
    power=np.pad(g.p_mid, [(int(g.nSlices/2)*n_pad, (g.nSlices-(int(g.nSlices/2))))*n_pad, (0, 0)], mode='constant')
    phase=np.pad(g.phi_mid, [(int(g.nSlices/2)*n_pad, (g.nSlices-(int(g.nSlices/2))))*n_pad, (0, 0)], mode='constant') #not supported by the numpy 1.6.2


    spectrum = abs(fft(np.sqrt( np.array(power)) * np.exp( 1.j* np.array(phase) ) , axis=0))**2/sqrt(g.nSlices)/(2*g.leng/g('ncar'))**2/1e10
    e_0=1239.8/g('xlamds')/1e9

    
    g.freq_ev1 = h * fftfreq(len(spectrum), d=g('zsep') * g('xlamds') / c)+e_0
    lamdscale=1239.8/g.freq_ev1

    lamdscale_array=np.swapaxes(np.tile(lamdscale,(g.nZ,1)),0,1)    
    
#    for std calculation
#    spectrum_lamdpos=np.sum(spectrum*lamdscale_array/np.sum(spectrum,axis=0),axis=0)
#    spectrum_lamdwidth=sqrt(np.sum(spectrum*(lamdscale_array-spectrum_lamdpos)**2/np.sum(spectrum,axis=0),axis=0))
    

    ax_spectrum.plot(fftshift(lamdscale), fftshift(spectrum[:,zi]), 'r-')
    ax_spectrum.text(0.5, 0.5,'on axis', horizontalalignment='center', verticalalignment='center')
    ax_spectrum.set_ylabel('P($\lambda$) [a.u.]')
    ax_spectrum.set_xlabel('$\lambda$ [nm]')
    ax_spectrum.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_spectrum.get_yaxis().get_major_formatter().set_scientific(True)
    ax_spectrum.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))#[:,75,75]
    if np.amin(lamdscale) != np.amax(lamdscale):
        ax_spectrum.set_xlim([np.amin(lamdscale), np.amax(lamdscale)])
    ax_phase.set_xlabel('s [$\mu$m]')


    maxspectrum_index=np.argmax(spectrum[:,zi])
    maxspectrum_wavelength=lamdscale[maxspectrum_index]*1e-9

    phase=unwrap(g.phi_mid[:,zi])
    
    phase_cor=np.arange(g.nSlices)*(maxspectrum_wavelength-g('xlamds'))/g('xlamds')*g('zsep')*2*pi
    phase_fixed=phase+phase_cor
    n=1
    phase_fixed = ( phase_fixed + n*pi) % (2 * n*pi ) - n*pi
    ax_phase.plot(s, phase_fixed, 'k-',linewidth=0.5)
    ax_phase.set_ylabel('$\phi$ [rad]')
    ax_phase.set_ylim([-pi, pi])

    number_ticks=6

    # ax_spectrum.yaxis.major.locator.set_params(nbins=number_ticks)
    
    ax_phase.xaxis.major.locator.set_params(nbins=number_ticks)
    ax_power.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_spectrum.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_bunching.yaxis.major.locator.set_params(nbins=number_ticks)
    ax_curr.yaxis.major.locator.set_params(nbins=number_ticks)

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

    ax_bunching.tick_params(axis='y', which='both', colors='grey')
    ax_bunching.yaxis.label.set_color('grey')  
    
    ax_power.tick_params(axis='y', which='both', colors='g')
    ax_power.yaxis.label.set_color('g')    
    ax_power.yaxis.get_offset_text().set_color(ax_power.yaxis.label.get_color())
    ax_spectrum.yaxis.get_offset_text().set_color(ax_spectrum.yaxis.label.get_color())
    

    if save!=False:
        if save==True:
            save='png'
        fig.savefig(g.path+'_z_'+str(z)+'m.'+str(save),format=save)
    
    return fig



def gen_outplot_scanned_z(g, figsize=(8, 10), legend = True, fig_name = None, z=inf, save=False):
#    max_yticks = 7
    if g('itdp')==True:
        print('    plotting scan at '+str(z)+' [m]')
        print('!     Not implemented yet for time dependent, skipping')
        return
        
    if g('iscan')==0 and g('scan')==0:
        print('    plotting scan at '+str(z)+' [m]')
        print('!     Not a scan, skipping')
        return
    
    import matplotlib.ticker as ticker
    
    if z==inf:
#        print 'Showing profile parameters at the end of undulator'
        z=np.amax(g.z)

    elif z>np.amax(g.z):
#        print 'Z parameter too large, setting to the undulator end'
        z=np.amax(g.z)
    
    elif z<np.amin(g.z):
#        print 'Z parameter too small, setting to the undulator entrance'    
        z=np.amin(g.z)
        
    zi=np.where(g.z>=z)[0][0]
    z=g.z[zi];
    
    print('    plotting scan at '+str(z)+' [m]')



    font_size = 1
    if fig_name is None:
        if g.filename is '':
            fig = plt.figure('Genesis scan at '+str(z)+'m')
        else:
            fig = plt.figure('Genesis scan at '+str(z)+'m '+g.filename)
    else:
        fig = plt.figure(fig_name)
    fig.set_size_inches(figsize,forward=True)
    plt.rc('axes', grid=True)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    # left, width = 0.1, 0.85
    
    plt.clf()
    
    ax_curr=fig.add_subplot(2, 1, 1)
    ax_curr.clear()
    ax_energy=fig.add_subplot(2, 1, 2,sharex=ax_curr)
    ax_energy.clear()    
#    ax_phase=fig.add_subplot(4, 1, 3,sharex=ax_curr)
#    ax_phase.clear()
#    ax_spectrum=fig.add_subplot(4, 1, 4)
#    ax_spectrum.clear()

    for ax in [ax_curr]:#, ax_energy: #ax_phase, ax_spectrum, 
        for label in ax.get_xticklabels():
            label.set_visible(False)

    #
    
    
    fig.subplots_adjust(hspace=0)
    
    s=g.scv #scan value is written to current colunm
        
    
    ax_curr.plot(s, np.linspace(g('curpeak'),g('curpeak'),len(s)), 'k--')
    ax_curr.set_ylabel('I[kA]')


    ax_power = ax_curr.twinx()
    ax_power.grid(False)
    ax_power.plot(s,g.p_int[:,zi],'g-',linewidth=1.5)    
    ax_power.set_ylabel('Power [W]')
    ax_power.set_ylim([0, np.amax(g.p_int[:,zi])])
    ax_power.get_yaxis().get_major_formatter().set_useOffset(False)
    ax_power.get_yaxis().get_major_formatter().set_scientific(True)
    ax_power.get_yaxis().get_major_formatter().set_powerlimits((-3, 4))#[:,75,75]
    
#    ax_power.get_xaxis().get_offset_text().set_x(1.1)

    ax_energy.plot(s, g.el_energy[:,zi]*0.511e-3, 'b-', s, (g.el_energy[:,zi]+g.el_e_spread[:,zi])*0.511e-3, 'r--',s, (g.el_energy[:,zi]-g.el_e_spread[:,zi])*0.511e-3, 'r--')
    ax_energy.set_ylabel('$E\pm\sigma_E$\n[GeV]')
#    ax_energy.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3), useOffset=False)
    ax_energy.ticklabel_format(useOffset=False, style='plain')   
    ax_energy.get_xaxis().get_major_formatter().set_useOffset(False)
    ax_energy.get_xaxis().get_major_formatter().set_scientific(True)


    ax_bunching = ax_energy.twinx()
    ax_bunching.plot(s,g.bunching[:,zi],'grey',linewidth=0.5)
    ax_bunching.set_ylabel('Bunching')
    ax_bunching.grid(False)
    
    
#    ax_power.yaxis.major.locator.set_params(nbins=number_ticks)
#    ax_energy.yaxis.major.locator.set_params(nbins=number_ticks)

    # ax_energy.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1e'))

    plt.xlim(s[0], s[-1])

    fig.subplots_adjust(top=0.95, bottom=0.2, right=0.85, left=0.15)

    #fig.set_size_inches((8,8),forward=True)
    
   
    ax_energy.tick_params(axis='y', which='both', colors='b')
    ax_energy.yaxis.label.set_color('b')    

    ax_bunching.tick_params(axis='y', which='both', colors='grey')
    ax_bunching.yaxis.label.set_color('grey')
    
    ax_power.tick_params(axis='y', which='both', colors='g')
    ax_power.yaxis.label.set_color('g')    
    ax_power.yaxis.get_offset_text().set_color(ax_power.yaxis.label.get_color())
    
    if save!=False:
        if save==True:
            save='png'
        fig.savefig(g.path+'_z_'+str(z)+'m_scan.'+str(save),format=save)
    
    return fig


def gen_outplot(handle=None,save='eps',show=False,debug=0):
    #picks as an input "GenesisOutput" object, file path of directory as strings. 
    #plots e-beam evolution, radiation evolution, initial and final simulation window
    #If folder path is provided, all *.gout and *.out files are plotted
    import os
    from ocelot.adaptors.genesis import readGenesisOutput, GenesisOutput
    
    plt.ioff()
    
    
    if os.path.isdir(str(handle)):
        handles=[]
        for root, dirs, files in os.walk(handle):
            for name in files:
                if name.endswith('.gout') or name.endswith('.out'):
                    handles.append(os.path.join(root, name))
        print('\n  plotting all files in '+str(handle))
    else:
        handles=[handle]
    
    for handle in handles:
    
        if os.path.isfile(str(handle)):
            handle=readGenesisOutput(handle,readall=1,debug=debug)
            
        if isinstance(handle,GenesisOutput):
            f1=gen_outplot_e(handle,save=save)
            f2=gen_outplot_ph(handle,save=save)
            f3=gen_outplot_z(handle, z=0,save=save)
            f4=gen_outplot_z(handle, z=inf,save=save)
    
    if show==True:
        print('    showing plots, close all to proceed')
        plt.show()
        
    if save!=False:
        print('    plots recorded to *.'+str(save)+' files')
    
    return [f1,f2,f3,f4]



def fwhm3(valuelist, peakpos=-1):
    """calculates the full width at half maximum (fwhm) of some curve.
    the function will return the fwhm with sub-pixel interpolation. It will start at the maximum position and 'walk' left and right until it approaches the half values.
    INPUT: 
    - valuelist: e.g. the list containing the temporal shape of a pulse 
    OPTIONAL INPUT: 
    -peakpos: position of the peak to examine (list index)
    the global maximum will be used if omitted.
    OUTPUT:
    -fwhm (value)
    """
    if peakpos== -1: #no peakpos given -> take maximum
        peak = np.max(valuelist)
        peakpos = np.min( np.nonzero( valuelist==peak  )  )

    peakvalue = valuelist[peakpos]
    phalf = peakvalue / 2.0

    # go left and right, starting from peakpos
    ind1 = peakpos
    ind2 = peakpos   

    while ind1>2 and valuelist[ind1]>phalf:
        ind1=ind1-1
    while ind2<len(valuelist)-1 and valuelist[ind2]>phalf:
        ind2=ind2+1  
    #ind1 and 2 are now just below phalf
    grad1 = valuelist[ind1+1]-valuelist[ind1]
    grad2 = valuelist[ind2]-valuelist[ind2-1]
    #calculate the linear interpolations
    p1interp= ind1 + (phalf -valuelist[ind1])/grad1
    p2interp= ind2 + (phalf -valuelist[ind2])/grad2
    #calculate the width
    width = p2interp-p1interp
    return (peakpos,width)

