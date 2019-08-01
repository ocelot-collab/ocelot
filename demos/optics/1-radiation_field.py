# -*- coding: utf-8 -*-
"""
Created on Mon Oct 8 2018
@author: sserkez
"""
import logging
import ocelot
from ocelot.common.globals import *
from ocelot.optics.wave import imitate_sase_dfl, wigner_dfl, dfl_waistscan, generate_gaussian_dfl
from ocelot.gui.dfl_plot import plot_dfl, plot_wigner, plot_dfl_waistscan
from copy import deepcopy


_logger = logging.getLogger(__name__)
#_logger.setLevel(logging.INFO)
_logger.setLevel(logging.DEBUG)


#%%

E_pohoton = 7000 #SVEA central photon energy [eV]

##simple field example
kwargs={'xlamds':(h_eV_s * speed_of_light / E_pohoton),
        'rho':1.0e-4, 
        'shape':(101,101,500), 
        'dgrid':(200e-6,200e-6,35e-6), 
        'power_rms':(10e-6,10e-6,3e-6), 
        'power_center':(0,0,None), 
        'power_angle':(0,0), 
        'power_waistpos':(0,0), 
        'wavelength':None, 
        'zsep':None, 
        'freq_chirp':0, #5e-3 
        'en_pulse':None, 
        'power':1e6,
        }


##complicated field example
#kwargs={'xlamds':(h_eV_s * speed_of_light / E_pohoton),
#        'rho':1.0e-4, 
#        'shape':(101,91,500), 
#        'dgrid':(200e-6,250e-6,35e-6), 
#        'power_rms':(10e-6,15e-6,3e-6),
#        'power_center':(15e-6,-30e-6,15e-6), 
#        'power_angle':(0,10e-6), 
#        'power_waistpos':(-2,-4), 
#        'wavelength':None, 
#        'zsep':None, 
#        'freq_chirp':1e-8, 
#        'en_pulse':None, 
#        'power':1e6,
#        }


        

##generating RadiationField() objects
#help(generate_dfl)
dfl = generate_gaussian_dfl(**kwargs);  #Gaussian
#dfl = imitate_sase_dfl(**kwargs); help(imitate_sase_dfl) #SASE-like (same as Gaussian, but modulated in time by amplitude and phase to model Gaussian random process)

#help(plot_dfl)
plot_dfl(dfl,
         domains=None,
         z_lim=[],
         xy_lim=[],
         figsize=4,
         cmap='viridis', #jet
         phase=False,
         fig_name='default_dfl_plot',
         auto_zoom=False,
         column_3d=True,
         savefig=False,
         showfig=True,
         return_proj=False,
         line_off_xy=True,
         log_scale=0,
         debug=1,
         cmap_cutoff=0)


plot_dfl(dfl, fig_name='dfl_xy', column_3d=0, phase=1)
plot_dfl(dfl, fig_name='dfl', domains='kt')
plot_dfl(dfl, fig_name='dfl', domains='sf')
plot_dfl(dfl, fig_name='dfl', domains='kf')
plot_dfl(dfl, fig_name='dfl', log_scale=1)

#calculating and plotting Wigner distribution
wig = wigner_dfl(dfl, 
                 method='mp', 
                 pad=1)

plot_wigner(wig, 
            x_units='fs', #fs
            y_units='ev', #nm
            x_lim=(None, None), 
            y_lim=(None, None), 
            downsample=1, 
            autoscale=None, 
            figsize=4, 
            cmap='seismic', 
            fig_name=None, 
            savefig=False, 
            showfig=True, 
            plot_proj=1, 
            plot_text=1, 
            plot_moments=0, #1
            )

#%% propagating

dfl = generate_gaussian_dfl(**kwargs)

#dfl.curve_wavefront(r=np.inf,
#                    plane='xy',
#                    domain_z=None,
#                    )

plot_dfl(dfl, fig_name='dfl_before_prop')

dfl0 = deepcopy(dfl)
dfl0.prop(z=5, #propagation distance [m]
          )
plot_dfl(dfl0, fig_name='dfl_prop')

dfl1 = deepcopy(dfl)
dfl1.prop_m(z=10,
            m=2, #magnification of the transverse grid size
            ) 
plot_dfl(dfl1, fig_name='dfl_prop_m')

#%% Scanning for waist position

w_scan = dfl_waistscan(dfl, 
                       z_pos=np.linspace(-20,20,15), 
                       projection=0)

plot_dfl_waistscan(w_scan, 
                   fig_name='waist scan results', 
                   showfig=True, 
                   savefig=False)

#%% propagation beyond paraxial approximation

kwargs={'xlamds':(h_eV_s * speed_of_light / E_pohoton),
        'rho':1.0e-4, 
        'shape':(101,101,50), 
        'dgrid':(100e-9,100e-9,10e-9), 
        'power_rms':(1e-9,1e-9,0.1e-9), 
        'power_center':(0,0,None), 
        'power_angle':(0,0), 
        'power_waistpos':(0,0), 
        'wavelength':None, 
        'zsep':None, 
        'freq_chirp':0, 
        'en_pulse':None, 
        'power':1e6,
        }

dfl = generate_gaussian_dfl(**kwargs);
plot_dfl(dfl, fig_name='dfl_non_paraxial')
#plot_dfl(dfl, domains='sf', fig_name='dfl_non_paraxial')
#plot_dfl(dfl, domains='kt', fig_name='dfl_non_paraxial')

dfl1 = deepcopy(dfl)
dfl1.prop_m(z=0.2,
            m=200000, #magnification of the transverse grid size
            fine=1 #if fine==1 every frequency component is propagated independently
            ) 
plot_dfl(dfl1, fig_name='dfl_non_paraxial_prop_spherical_front')