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
params = {'backend': 'ps', 'axes.labelsize': 15, 'font.size': 15, 'legend.fontsize': 24, 'xtick.labelsize': 19,  'ytick.labelsize': 19, 'text.usetex': True}
rcParams.update(params)
rc('text', usetex=True) # required to have greek fonts on redhat

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 15}

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
    ax_spectrum.text(0, 0,r"on axis", fontsize=15)#horizontalalignment='center', verticalalignment='center',
    ax_spectrum.set_ylabel('P($\lambda$) [a.u.]')
    ax_spectrum.set_xlabel('$\lambda$ [nm]')
    ax_spectrum.set_ylim(ymin=0)
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


def gen_outplot(handle=None,save='png',show=False,debug=0):
    #picks as an input "GenesisOutput" object, file path of directory as strings. 
    #plots e-beam evolution, radiation evolution, initial and final simulation window
    #If folder path is provided, all *.gout and *.out files are plotted
    import os
    from ocelot.adaptors.genesis import GenesisOutput, readGenesisOutput, readRadiationFile
    
    plt.ioff()
    
    if save==True:
        save='png'
    
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
        
        if os.path.isfile(handle.path+'.dfl'):
            dfl=readRadiationFile(handle.path+'.dfl', handle.ncar)
            f5=plot_dfl(dfl, handle,save=save)
            
    if show==True:
        print('    showing plots, close all to proceed')
        plt.show()
        
    if save!=False:
        print('    plots recorded to *.'+str(save)+' files')
    
    # return [f1,f2,f3,f4]


def plot_dfl(dfl, g, figsize=3, legend = True, fig_name = None, save=False):
    
    print('    plotting dfl file')
    
    if dfl.shape[0]!=1:
        column_3d=True
    else:
        column_3d=False
    
    dfl=swapaxes(dfl[::-1,:,:],2,1) # zyx -> zxy
    
    #number of mesh points
    ncar_x=dfl.shape[1]
    leng_x=g.leng*1e6 #transverse size of mesh [um], to be upgraded
    ncar_y=dfl.shape[2]
    leng_y=g.leng*1e6
    if column_3d==True:
        ncar_z=dfl.shape[0]
        leng_z=(max(g.s)-min(g.s))*1e6
    
    x = np.linspace(-leng_x/2, leng_x/2, ncar_x)
    y = np.linspace(-leng_y/2, leng_y/2, ncar_y)
    z = np.linspace(0, leng_z, ncar_z)
    
    if fig_name is None:
        if g.filename is '':
            fig = plt.figure('Radiation field')
        else:
            fig = plt.figure('Radiation field, '+g.filename)
    else:
        fig = plt.figure(fig_name)
    fig.clf()
    fig.set_size_inches(((3+2*column_3d)*figsize,3*figsize),forward=True)
    # plt.rc('axes', grid=True)
    # plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    
    cmap_int = plt.get_cmap('jet') #change to convenient
    cmap_ph = plt.get_cmap('hsv')
    
    #calculate transverse projection, remove z dimention
    
    dfl_int=abs(dfl)**2
    xy_proj=sqrt((dfl_int).sum(0))*exp(1j*angle(dfl.sum(0))) #(amplitude-like) view from front, sum of square of amplitudes with phase as sum of phasors (latter is dedicated for illustration purposes: good to see an averaged wavefront)
    # xy_proj=sum(dfl,0); #view from front
    yz_proj=sum(dfl_int,1); #intensity view from side
    xz_proj=sum(dfl_int,2); #intensity view from top
    del dfl_int

    
    # x_line=xy_proj[]
    # y_line=xy_proj[]
    
    int_proj=abs(xy_proj)**2
    ph_proj=angle(xy_proj)
    
    x_proj=sum(int_proj,1)
    y_proj=sum(int_proj,0)
    
    x_line=int_proj[:,int((ncar_y-1)/2)]
    y_line=int_proj[int((ncar_x-1)/2),:]
    
    if max(x_line)!=0 and max(y_line)!=0:
        x_line,y_line=x_line/max(x_line),y_line/max(y_line)
    
    
    
    #X=sqrt(sum(abs(X).^2,3)).*exp(1i.*angle(mean(X,3))); #%Matlab 2D field calculation
    # normI = BoundaryNorm(levelsI, ncolors=cmapI.N, clip=True)
    # normP = BoundaryNorm(levelsP, ncolors=cmapP.N, clip=True)
    
    
    
    
    
    
    ax_int=fig.add_subplot(2, 2+column_3d, 1)
    # ax_int.pcolormesh(x, y, int_proj, cmap=cmap_int)
    intplt=ax_int.pcolormesh(x, y, swapaxes(int_proj,1,0))
    ax_int.set_title('Intensity', fontsize=15)
    ax_int.axis('equal')
    # ax_int.axes.get_xaxis().set_visible(False)
    ax_int.set_xlabel(r'x [$\mu m$]')
    ax_int.set_ylabel(r'y [$\mu m$]')
    
    ax_ph=fig.add_subplot(2, 2+column_3d, 4+column_3d, sharex=ax_int,sharey=ax_int)
    # ax_ph.pcolormesh(x, y, ph_proj, cmap=cmap_ph)
    ax_ph.pcolormesh(x, y, swapaxes(ph_proj,1,0))
    #ax_ph.axis('equal')
    ax_ph.axis([min(x),max(x),min(y),max(y)])
    ax_ph.set_title('Phase', fontsize=15)
    # ax_ph.set_xlabel(r'[$\mu m$]')
    # ax_ph.set_ylabel(r'[$\mu m$]')
    
    ax_proj_x=fig.add_subplot(2, 2+column_3d, 3+column_3d, sharex=ax_int)
    ax_proj_x.plot(x,x_line)
    ax_proj_x.set_title('X projection', fontsize=15)
    x_line_f, fwhm_x=fwhm_gauss_fit(x,x_line)
    ax_proj_x.plot(x,x_line_f)
    ax_proj_x.text(0.95, 0.95,'FWHM= '+str(round_sig(fwhm_x,3))+r'$\mu m$', horizontalalignment='right', verticalalignment='top', transform = ax_proj_x.transAxes, fontsize=12)
    ax_proj_x.set_ylim(ymin=0,ymax=1)
    
    
    ax_proj_y=fig.add_subplot(2, 2+column_3d, 2, sharey=ax_int)
    ax_proj_y.plot(y_line,y)
    ax_proj_y.set_title('Y projection', fontsize=15)
    y_line_f, fwhm_y=fwhm_gauss_fit(y,y_line)
    ax_proj_y.plot(y_line_f,y)
    ax_proj_y.text(0.95, 0.95,'FWHM= '+str(round_sig(fwhm_y,3))+r'$\mu m$', horizontalalignment='right', verticalalignment='top', transform = ax_proj_y.transAxes, fontsize=12)
    ax_proj_y.set_xlim(xmin=0,xmax=1)

    
    
    if column_3d:
        ax_proj_xz=fig.add_subplot(2, 2+column_3d, 6)
        ax_proj_xz.pcolormesh(z, x, swapaxes(xz_proj,1,0))
        ax_proj_xz.set_title('top view', fontsize=15)
        ax_proj_xz.axis('tight')
        ax_proj_xz.set_xlabel(r'z [$\mu m$]')
        ax_proj_yz=fig.add_subplot(2, 2+column_3d, 3,sharey=ax_int,sharex=ax_proj_xz)
        ax_proj_yz.pcolormesh(z, y, swapaxes(yz_proj,1,0))
        ax_proj_yz.set_title('side view', fontsize=15)
        ax_proj_yz.axis('tight')

        
        
        
    cbar=0
    if cbar:
        fig.subplots_adjust(top=0.95, bottom=0.05, right=0.85, left=0.1)
        #fig.subplots_adjust()
        cbar_int = fig.add_axes([0.89, 0.15, 0.015, 0.7])
        cbar=plt.colorbar(intplt, cax=cbar_int)# pad = -0.05 ,fraction=0.01)
        # cbar.set_label(r'[$ph/cm^2$]',size=10)
        cbar.set_label(r'a.u.',size=10)
    
    
    
    
    # ax_int.get_yaxis().get_major_formatter().set_useOffset(False)
    # ax_int.get_yaxis().get_major_formatter().set_scientific(True)
    # ax_ph.get_yaxis().get_major_formatter().set_useOffset(False)
    # ax_ph.get_yaxis().get_major_formatter().set_scientific(True)
    
    ax_int.set_aspect('equal')
    ax_int.autoscale(tight=True)
    
    subplots_adjust(wspace=0.4,hspace=0.4)
    
    
    if save!=False:
        if save==True:
            save='png'
        fig.savefig(g.path+'_dfl.'+str(save),format=save)
        
    return fig
    
    

# def plot_dfl_o(g, figsize=(8, 10), legend = True, fig_name = None, save=False):
# # def plotfield1(name ,Xs, Xf, nx, Ys, Yf, ny, arI, arP, E_ph):
    # import matplotlib.pyplot as plt
    # from matplotlib.colors import BoundaryNorm
    # from matplotlib.ticker import MaxNLocator
    # import numpy as np
    # import math
    # import pylab as pl
    # from pylab import clf, axis

    # # generate 2 2d grids for the x & y bounds

    # x = 1e6 * np.linspace(Xs, Xf, nx)
    # y = 1e6 * np.linspace(Ys, Yf, ny)

    # dataI = np.array( arI )/E_ph/1.6e-19/10000
    # #dataI =x/E_ph/1.6e-19 for x in dataI
    # zI=dataI.reshape( ny, nx )

    # dataP = np.array( arP )
    # zP=dataP.reshape( ny, nx )

    # #x and y are bounds, so z should be the value *inside* those bounds.
    # # Therefore, remove the last value from the z array.
    # #zI = zI[:-1, :-1]
    # levelsI = MaxNLocator(nbins=50).tick_values(zI.min(), zI.max())
    # #zP = zP[:-1, :-1]
    # levelsP = MaxNLocator(nbins=50).tick_values(-math.pi, math.pi)


    # # pick the desired colormap, sensible levels, and define a normalization
    # # instance which takes data values and translates those into levels.
    # cmapI = plt.get_cmap('jet')
    # cmapP = plt.get_cmap('hsv')
    # normI = BoundaryNorm(levelsI, ncolors=cmapI.N, clip=True)
    # normP = BoundaryNorm(levelsP, ncolors=cmapP.N, clip=True)

    # fig = plt.figure(name)
    # clf()

    # # plt.subplot(2, 2, 1)-------------------------------------
    # p1=fig.add_subplot(2, 2, 1)
    # p1.clear()

    # pp1=p1.pcolormesh(x, y, zI, cmap=cmapI)#, norm=normI)

    # #plt.colorbar(pp1,shrink=0.9,fraction=0.01)# pad = -0.05)
    # # # set the limits of the plot to the limits of the data
    # p1.axis([x.min(), x.max(), y.min(), y.max()])
    # #p1.title('intensity')
    # plt.title('Intensity')
    # axis('equal')

    # p1.set_xlabel(r'[$\mu m$]')
    # p1.set_ylabel(r'[$\mu m$]')
    # # axis([Xs*1e6, Xf*1e6, Ys*1e6, Yf*1e6])

    # # plt.subplot(2, 2, 2)--------------------------------------
    # p2=fig.add_subplot(2, 2, 2, sharey=p1)
    # p2.clear()

    # #arIm=reshape(np.matrix(arI), [nx, ny])
    # #arIpx=np.sum(arIm,1)
    # #plt.plot(np.linspace(wfr.mesh.xStart,wfr.mesh.xFin,wfr.mesh.nx),arIpx)

    # arIpy=np.sum(zI,1)
    # arIpy=arIpy/max(arIpy)
    # #
    # #p2.plot(x, arIpx)
    # p2.plot(arIpy,y)
    # arIpyf, fwhmy=fwhm_gauss_fit(y,arIpy)
    # p2.plot(arIpyf,y)

    # pl.text(0.95, 0.95,'FWHM= '+str(round_sig(fwhmy,3))+r'$\mu m$',
    # horizontalalignment='right',
    # verticalalignment='top',
    # transform = p2.transAxes)  

    # # plt.subplot(2, 2, 3)--------------------------------------
    # p3=fig.add_subplot(2, 2, 3, sharex=p1)
    # p3.clear()

    # #arIm=reshape(np.matrix(arI), [nx, ny])
    # #arIpx=np.sum(arIm,1)
    # #plt.plot(np.linspace(wfr.mesh.xStart,wfr.mesh.xFin,wfr.mesh.nx),arIpx)

    # arIpx=np.sum(zI,0)
    # arIpx=arIpx/max(arIpx)

    # #p2.plot(x, arIpx)
    # p3.plot(x,arIpx)
    # arIpxf, fwhmx=fwhm_gauss_fit(x,arIpx)
    # #fwhmv=fwhm_gauss_fit(x,arIpx)[1]
    # #np.disp(fwhmv)
    # p3.plot(x,arIpxf)

    # pl.text(0.95, 0.95,'FWHM= '+str(round_sig(fwhmx,3))+r'$\mu m$',
    # horizontalalignment='right',
    # verticalalignment='top',
    # transform = p3.transAxes)

    # # # set the limits of the plot to the limits of the data
    # #p2.axis([x.min(), x.max()])
    # # plt.title('intensity')
    # #axis('equal')
    # # axis([Xs*1e6, Xf*1e6, Ys*1e6, Yf*1e6])

    # # plt.subplot(2, 2, 4)---------------------------------------
    # p4=fig.add_subplot(2, 2, 4, sharex=p1)
    # p4.clear()

    # p4.pcolormesh(x, y, zP, cmap=cmapP, norm=normP)
    # #plt.colorbar()
    # # set the limits of the plot to the limits of the data
    # # plt.axis([x.min(), x.max(), y.min(), y.max()])
    # plt.title('Phase')
    # axis('equal')
    # # axis([Xs*1e6, Xf*1e6, Ys*1e6, Yf*1e6])

    # ## plt.subplot(2, 2, 4)---------------------------------------
    # #pha=fig.add_subplot(2, 2, 4)
    # #pha.clear()
    # #
    # #plt.pcolormesh(x, y, zP, cmap=cmapP, norm=normP)
    # #plt.colorbar()
    # ## set the limits of the plot to the limits of the data
    # ## plt.axis([x.min(), x.max(), y.min(), y.max()])
    # #plt.title('Phase')
    # #axis('equal')
    # ## axis([Xs*1e6, Xf*1e6, Ys*1e6, Yf*1e6])

    # fig.subplots_adjust(top=0.95, bottom=0.05, right=0.85, left=0.1)
    # #fig.subplots_adjust()
    # cbar_int = fig.add_axes([0.89, 0.15, 0.015, 0.7])
    # cbar=plt.colorbar(pp1, cax=cbar_int)# pad = -0.05 ,fraction=0.01)
    # cbar.set_label(r'[$ph/cm^2$]',size=10)
    # plt.show()
    # fig.canvas.draw()
    # draw()
    # return fig


def round_sig(x, sig=2):
    from math import log10, floor
    return round(x, sig-int(floor(log10(x)))-1)
    
    
def fwhm_gauss_fit(X,Y):
    import numpy as np
    import scipy.optimize as opt
    
    def gauss(x, p): # p[0]==mean, p[1]==stdev p[2]==peak
        return p[2]/(p[1]*np.sqrt(2*np.pi))*np.exp(-(x-p[0])**2/(2*p[1]**2))
        
    p0 = [0,max(X)/2,max(Y)]
    errfunc = lambda p, x, y: gauss(x, p) - y
    p1, success = opt.leastsq(errfunc, p0[:], args=(X, Y))
    fit_mu,fit_stdev,ampl = p1
    Y1=gauss(X,p1)
    FWHM = 2*np.sqrt(2*np.log(2))*fit_stdev
    return (Y1, FWHM)

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

