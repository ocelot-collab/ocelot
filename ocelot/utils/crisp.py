# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 15:11:54 2019

Module for the different tasks that need to be done to reconstruct the current
profile from a formfactor measurement

@author: lockmann

Example:
---
from ocelot import *
import ocelot.utils.crisp as srisp
import matplotlib.pyplot as plt

sigma_tau = 10e-6
charge = 250e-12
parray_init = generate_parray(sigma_tau=sigma_tau, sigma_p=2.5 / 14000, chirp=-0.001, charge=charge,
                                          nparticles=200000, energy=14, tws=None, shape="tri")
s, I = get_current(parray_init)
# generate CRISP spectrum

frequencies, formfactors, formfactor_noise, detlim = srisp.get_crisp_signal(s, I, n_shots = 20, which_set="both")


recon_time, recon_current, t_rms = srisp.master_recon(frequencies, formfactors, formfactor_noise, detlim,
                                                charge=parray_init.total_charge, method='KKstart',
                                                        channels_to_remove = [], show_plots = False)
plt.figure(figsize=(12,5))
ax1 = plt.subplot(121)
plt.title('Original Current')
ax1.plot(s/3e8*1e15, I*1e-3, label = 'Input')
ax1.set_ylabel("I [kA]")
ax1.set_xlabel("s [fs]")

ax2 = plt.subplot(122)
plt.title('Spectrometer Formfactor')
ax2.errorbar(frequencies*1e-12, formfactors, yerr=formfactor_noise, fmt='o', label='|F|')
ax2.fill_between(frequencies*1e-12, np.zeros(np.size(detlim)), y2 = detlim, label = 'Noise Floor', alpha = 0.3, color='gray')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('f (THz)')
ax2.set_ylabel('|F|')
plt.legend()

plt.figure(40)
plt.title('Reconstruction')
plt.plot(recon_time*1e15, recon_current*1e-3, label='Reconstruction')
plt.plot(s/3e8*1e15, I*1e-3, label='Input')
plt.xlabel('t (fs)')
plt.ylabel('I (kA)')
plt.xlim([-5*t_rms*1e15, 5*t_rms*1e15])
plt.legend()
plt.show()

"""

import numpy as np
import matplotlib.pylab as plt
from scipy.io import loadmat
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.signal import savgol_filter

"""
General Helper Functions
"""


def first_moment(intensity, x):
    return np.sum(intensity * x) / np.sum(intensity)


def second_moment(intensity, x):
    return np.sum(intensity * (x - first_moment(intensity, x)) ** 2) / np.sum(intensity)


def super_gauss(x, amp, order, sigma):
    """
    Returns a super gaussian function with the given parameters
    """
    return amp * np.exp(-0.5 * (x / sigma) ** (2 * int(order)))


def gauss(x, sigma):
    """
    Gauss Function with amplitude 1
    """
    return super_gauss(x, 1, 1, sigma)


def gauss_wamp(x, amp, sigma):
    """
    Gauss Function with variable amplitude bigger than 1
    """
    if amp < 1:
        amp = 1
    return amp * super_gauss(x, 1, 1, sigma)


def nils_fft(xin, datain):
    """
    Function to return
    """
    fouriertransform = np.fft.fft(datain)
    quatschfaktor = 1 / np.sqrt(2. * np.pi) * (xin[1] - xin[0])
    fouriertransform = quatschfaktor * fouriertransform
    fx = np.fft.fftfreq(np.size(xin), d=xin[1] - xin[0])
    #fx = sort(fx)
    wx = 2 * np.pi * fx
    return wx, fouriertransform


"""
Functions for processing of CRISP data
"""


def remove_certain_channels(channel_list, dataset):
    """
    Remove certain channels from the formfactor (broken ones for example)
    """
    new_data = np.copy(dataset)
    new_data = np.delete(new_data, channel_list, axis=1)
    return new_data


def sort_data(frequencies, data_to_sort):
    """
    Sort data according to frequency. Very important for everything to come.
    """
    if np.size(np.shape(data_to_sort)) > 1:
        result = np.empty_like(data_to_sort)
        for i in np.arange(np.shape(data_to_sort)[0]):
            this_sorted = data_to_sort[i][np.argsort(frequencies)]
            result[i] = this_sorted
    return np.sort(frequencies), result


def remove_loners(sorted_frequencies, sorted_formfactors, sorted_noise, sorted_det_limit):
    """
    we don't want single channels to pop up in regions where there is not much signal.
    Therefore, we demand 4 neighbours above the significanc level and must not be nan
    """
    my_maske = np.ones(np.size(sorted_frequencies), dtype='bool')
    for channel in np.arange(np.size(sorted_formfactors)):
        this_values = sorted_formfactors[channel - 2:channel + 2]
        this_det_limits = sorted_det_limit[channel - 2:channel + 2]
        if np.any(this_values < this_det_limits) or np.any(np.isnan(this_values)):
            my_maske[channel] = False
        #edit 21.03.2020: alle unter dem detection limit, sind auch kleiner als das detection limit aber mit unendlich großem fehler
    new_sorted_formfactors = np.copy(sorted_formfactors)
    new_sorted_noise = np.copy(sorted_noise)
    new_sorted_formfactors[~my_maske] = 0.5 * sorted_det_limit[~my_maske]
    new_sorted_noise[~my_maske] = 0.5 * sorted_det_limit[~my_maske]
    return sorted_frequencies[my_maske], sorted_formfactors[my_maske], sorted_noise[my_maske]
    #return sorted_frequencies, new_sorted_formfactors, new_sorted_noise


#def remove_outliers(sorted_formfactors, n_neighbours = 3):
#    """
#    Check for unphysical jumps, if the values change more than 0.5/2 in the
#    neighbours range, the data point is removed
#    """
#    channels = np.arange(np.size(sorted_formfactors))[::-1]
#    changes= sorted_formfactors[1:]/sorted_formfactors[:-1]
#    low_jumps = channels[changes<0.5]
#    high_jumps = channels[changes>2.]
#    distances = n

def extrapolation_low(new_frequencies, sorted_frequencies, sorted_formfactors, sorted_formfactor_noise):
    """
    Extrapolation to low Frequencies. As soon as 3 channels have a formfactor bigger than ff_value,
    A Gaussian function will be fitted and this data replaced by the fit.
    Question: Allow Formfactors bigger than 1? -> No
    """

    #find the value above which it will be extrapolated
    ff_value = 0.2
    while ff_value < 0.9 and np.size(sorted_formfactors[sorted_formfactors > ff_value]) > 8:
        ff_value = ff_value + 0.1

    channels = np.arange(np.size(sorted_frequencies))
    #find where the first time 3 ff are smaller than ff_value
    smallff_channels = channels[sorted_formfactors < ff_value]
    #scan to fnd the first time 3 smaller 0.8
    n_hits = 0
    i = 1
    while n_hits < 2:
        if smallff_channels[i] == smallff_channels[i - 1] + 1:
            #print(smallff_channels[i])
            n_hits = n_hits + 1
        else:
            n_hits = 0
        #print(n_hits)
        i = i + 1
    n_channels_tmp = smallff_channels[i]
    # Long Bunches:
    if n_channels_tmp < 4:
        n_channels = 60
    else:
        n_channels = n_channels_tmp
    #print('stopped at ' + str(sorted_frequencies[n_channels]*1e-12))
    #n_channels = np.amax(channels[sorted_formfactors>0.8])
    fit_freqs, fit_formfactors = sorted_frequencies[:n_channels], sorted_formfactors[:n_channels]
    guess = 1e12
    popt, pcov = curve_fit(gauss, fit_freqs, fit_formfactors, p0=guess)
    #guess_width = np.sqrt(-sorted_frequencies[n_channels]**2/ 2 / np.log(sorted_formfactors[n_channels]))
    #guess = np.array([1,guess_width])
    #popt, pcov = curve_fit(gauss_wamp, fit_freqs, fit_formfactors, p0 = guess)
    #for long bunches to only extrapolate but not replace:
    if n_channels_tmp < 4:
        n_channels = 0
    added_freqs = new_frequencies[new_frequencies < sorted_frequencies[n_channels]]
    added_formfactor = gauss(added_freqs, popt)  #[0], popt[1])
    new_formfactors = np.append(added_formfactor, sorted_formfactors[n_channels:])
    new_noise = np.append(np.zeros_like(added_formfactor), sorted_formfactor_noise[n_channels:])
    new_freqs = np.append(added_freqs, sorted_frequencies[n_channels:])
    return new_freqs, new_formfactors, new_noise


def extrapolation_high(new_frequencies, sorted_frequencies, sorted_formfactors, sorted_formfactor_noise,
                       index_to_cut=1):
    """
    Extrapolation to high Frequencies. A Gaussian will be fitted to the value
    at the data point at the highest available frequency.
    """
    last_freq = sorted_frequencies[-index_to_cut]
    last_mess = sorted_formfactors[-index_to_cut]
    high_sig_f = np.sqrt(-last_freq ** 2 / 2 / np.log(last_mess))

    fit_freqs = new_frequencies[new_frequencies > last_freq]
    fit_formfactors = gauss(fit_freqs, high_sig_f)
    new_formfactors = np.append(sorted_formfactors, fit_formfactors)
    new_noise = np.append(sorted_formfactor_noise, np.zeros_like(fit_formfactors))
    new_freqs = np.append(sorted_frequencies, fit_freqs)
    return new_freqs, new_formfactors, new_noise


def smooth_formfactors(formfactors, window_size=9):
    if window_size > 4:
        polyorder = 3
    elif window_size == 3:
        polyorder = 2
    else:
        polyorder = 0
    new_formfactors = savgol_filter(formfactors, window_size, polyorder)
    return new_formfactors


def make_linear_spacing(new_frequencies, sorted_frequencies, sorted_formfactors, sorted_formfactors_noise):
    #Interpolation
    f = interpolate.interp1d(sorted_frequencies, sorted_formfactors, bounds_error=True)
    new_formfactors = f(new_frequencies)
    g = interpolate.interp1d(sorted_frequencies, sorted_formfactors_noise, bounds_error=True)
    new_noise = g(new_frequencies)
    return new_formfactors, new_noise


"""
Models to fit the formfactors to. 
"""


def rechteck(x, sigma):
    data = super_gauss(x, 1, 50, sigma)
    return data


def unsym_gauss(x, sigma1, sigma2):
    data = np.zeros(np.size(x))
    data[x >= 0] = gauss(x[x >= 0], sigma1)
    data[x < 0] = gauss(x[x < 0], sigma2)
    return data


def cauchy(x, sigma):
    return 1 / (1 + (x / sigma) ** 2)


def unsym_cauchy(x, sigma1, sigma2):
    data = np.zeros(np.size(x))
    data[x >= 0] = cauchy(x[x >= 0], sigma1)
    data[x < 0] = cauchy(x[x < 0], sigma2)
    return data


def saw_tooth(x, sigma):
    data = 1 - 0.5 / sigma * x
    data[x < 0] = 0
    data[data < 0] = 0
    return (data)


def triangular(x, sigma1, sigma2):  #
    data = np.zeros(np.size(x))
    data[x >= 0] = 1 - 0.5 / sigma1 * x[x >= 0]
    data[x < 0] = 1 + 0.5 / sigma2 * x[x < 0]
    data[data < 0] = 0
    return (data)


def exp_step(x, sigma_t):
    expo = np.exp(-x / sigma_t)
    expo[x < 0] = 0
    return expo


#List of the models with proper start paras
guess_width = 100e-15
modellist = [[gauss, [guess_width]],
             [rechteck, [guess_width]],
             [unsym_gauss, [guess_width, guess_width]],
             [unsym_cauchy, [guess_width, guess_width]],
             [saw_tooth, [guess_width]],
             [triangular, [guess_width, guess_width]],
             [exp_step, [guess_width]],
             ]

"""
To find good model start
"""


def ff_from_model(model, freqs, modelparas):
    """
    Returns the normalised imaginary formfactor of the model
    """
    #print( model )
    global wx
    zeit = np.fft.fftfreq(2 * np.size(freqs) - 1, d=freqs[1] - freqs[0])
    zeit = np.sort(zeit)
    if np.size(modelparas) == 2:
        zeitdom = model(zeit, modelparas[0], modelparas[1])
    else:
        zeitdom = model(zeit, modelparas)
    #normalize to densitiy 1
    zeitdom = zeitdom / np.trapz(zeitdom, x=zeit)
    #print (modelparas)
    #zeitdom= model(zeit, *[modelparas])
    #ff_model = np.fft.fft(zeitdom)
    wx, ff_model = nils_fft(zeit, zeitdom)
    ff_model = ff_model[:np.size(freqs)] * np.sqrt(2 * np.pi)
    #ff_abs_model = np.abs(ff_model)
    #ff_model = ff_model / ff_abs_model[0]
    return ff_model


def fit_model(modeldata, freqs, formfactors, highest_freq):
    """
    Function to fit a model such that the formfactor modulus fits best to the
    data up to highest freq
    """
    #global model, model_guess
    model = modeldata[0]
    model_guess = modeldata[1]
    if np.size(model_guess) == 2:
        def this_function_to_fit(freqs, sigma1, sigma2):
            sigmas = [sigma1, sigma2]
            this_ff = ff_from_model(model, freqs, sigmas)
            return np.abs(this_ff)
    else:
        def this_function_to_fit(freqs, sigma):
            this_ff = ff_from_model(model, freqs, sigma)
            return np.abs(this_ff)

    formfactors_to_fit = formfactors[freqs < highest_freq]
    freqs_to_fit = freqs[freqs < highest_freq]

    popt, pcov = curve_fit(this_function_to_fit, freqs_to_fit, formfactors_to_fit, p0=model_guess)

    #Get goodness of fit
    #if np.size(popt) ==2:
    #    resulting_ff_abs = np.abs(ff_from_model(model, freqs_to_fit, popt[0], popt[1]))
    #else:
    resulting_ff_abs = np.abs(ff_from_model(model, freqs_to_fit, popt))
    criterium = np.sum((resulting_ff_abs - formfactors_to_fit) ** 2)
    return criterium, popt


def find_best_modelff(modellist, freqs, formfactors, highest_freq):
    """
    Find the best model from the list of models and returns that formfactor
    """
    n_models = np.shape(modellist)[0]
    gueten = []
    weiten = []
    for n in np.arange(n_models):
        this_res = fit_model(modellist[n], freqs, formfactors, highest_freq)
        gueten = np.append(gueten, this_res[0])
        weiten.append(this_res[1])
    #get the best suited one
    best_model = np.argmin(gueten)
    print('Best found model is ' + modellist[best_model][0].__name__)

    #get that models formfactor and return it
    return ff_from_model(modellist[best_model][0], freqs, weiten[best_model])


"""
GERCHBERG-SAXTON RECONSTRUCTION ALGORITHM
"""


def gerchberg_saxton(formfactors, formfactors_noise, start_phase, phase_condition=np.pi / 0.01):
    """
    Goes through the iterative phase adjustment until the phase differene is
    below phase_condition
    """
    j = complex(0, 1)
    seed_phase = start_phase
    difference = 1e6
    new_formfactor = np.copy(formfactors)
    n = 0
    nmax = 20
    while n < nmax:
        #while difference > phase_condition:
        prev_phase = np.copy(seed_phase)
        #Nur Formfaktoren austauschen, die nicht im Rauschen liegen
        new_formfactor = np.abs(new_formfactor)
        bools_too_low = new_formfactor < (formfactors - formfactors_noise)
        bools_too_high = new_formfactor > (formfactors + formfactors_noise)
        new_formfactor[bools_too_low] = formfactors[bools_too_low]
        new_formfactor[bools_too_high] = formfactors[bools_too_high]

        comp_formfactor = new_formfactor * np.exp(j * seed_phase)
        #erweitern, damit auch negative frequenzen etc. dabei sind. Diese Sortierung ist
        # wichtig, damit ein reales Zeitprofil rauskommt
        comp_formfactor = np.append(comp_formfactor, np.conj(comp_formfactor[1::])[::-1])

        #Fourier Transform to get current profile
        this_curr_prof = np.fft.ifft(comp_formfactor)
        this_curr_prof = np.real(this_curr_prof)
        this_curr_prof = np.roll(this_curr_prof, - int(np.argmax(this_curr_prof) - np.size(this_curr_prof) / 2))
        #plt.figure('Iteration')
        #plt.plot(this_curr_prof, label = str(n))
        #plt.legend()
        begin_curr_prof = np.copy(this_curr_prof)
        #Error-Reduction-Method: Set to 0 where violating
        #a) negative charges
        #violation = np.any(this_curr_prof<-1e-2)
        this_curr_prof[this_curr_prof < 0] = 0
        #b) being larger than 1ps
        #center = first_moment(this_curr_prof, zeit)
        #this_curr_prof[:100] = 0
        #this_curr_prof[-100:] = 0
        #c) zusammenhängend
        # get rid of postoscillations...
        #check low
        x = np.arange(np.size(this_curr_prof) / 2)
        low_positions = x[this_curr_prof[:np.size(x)] == 0]
        high_positions = x[this_curr_prof[np.size(x) - 1:] == 0]
        if np.size(low_positions) == 0:
            low_index = 0
        else:
            low_index = low_positions[-1]
        if np.size(high_positions) == 0:
            high_index = x[-1]
        else:
            high_index = high_positions[0] + x[-1]
        this_curr_prof[:int(low_index)] = 0
        this_curr_prof[int(high_index):] = 0
        #..done

        #plt.plot(this_curr_prof, label=step)
        #FFT to get new phase in time domain
        x = np.arange(np.size(this_curr_prof))
        """
        Change here
        """
        this_curr_prof = this_curr_prof  #/np.trapz(this_curr_prof, x = x)
        """
        """
        #new_formfactor = np.fft.fft(this_curr_prof)
        wx, new_formfactor = nils_fft(x, this_curr_prof)
        new_formfactor = new_formfactor * np.sqrt(2 * np.pi)
        new_formfactor = new_formfactor[:int((np.size(new_formfactor) + 1) / 2)]
        seed_phase = np.angle(new_formfactor)
        difference = np.sum(np.abs(seed_phase - prev_phase) ** 2)
        #difference = np.sum(np.abs(np.abs(new_formfactor[:20])-np.abs(inter_formfactors)))
        #Wenn ich so aufhöre, stimmt der Formfakor bei hohen Frequenzen nicht mehr überein:
        #if difference <=phase_condition:
        if n == nmax - 1:
            clean_curr_prof = np.copy(this_curr_prof)  # mit Cleanup
            this_curr_prof = begin_curr_prof  # ohne Clean up
            x = np.arange(np.size(this_curr_prof))
            #outcommented this
            #this_curr_prof = this_curr_prof /np.trapz(this_curr_prof, x = x)
            #new_formfactor = np.fft.fft(this_curr_prof)
            #wx, new_formfactor = nils_fft(x,this_curr_prof)
            #new_formfactor = new_formfactor * np.sqrt(2*np.pi)
            #new_formfactor = new_formfactor[:int((np.size(new_formfactor)+1)/2)]
        n = n + 1

    #this_curr_prof = np.roll(this_curr_prof, - int(np.argmax(this_curr_prof)-np.size(this_curr_prof)/2))
    #center fine
    #center = first_moment(this_curr_prof, np.arange(np.size(this_curr_prof)))
    center = first_moment(clean_curr_prof, np.arange(np.size(clean_curr_prof)))
    final_curr_prof = np.roll(this_curr_prof, -int(center - np.size(this_curr_prof) / 2))
    clean_curr_prof = np.roll(clean_curr_prof, -int(center - np.size(clean_curr_prof) / 2))
    return clean_curr_prof  #, final_curr_prof # # mit und ohne clean up


def kramers_kronig_phase(frequencies, mess_formfactors):
    """
    calculate the kramers kronig phase. Input frequencies and formfactors must be sorted
    """
    #Ln works only for positiv Formfactors
    these_freqs = frequencies[mess_formfactors > 0]
    these_ff = mess_formfactors[mess_formfactors > 0]
    indices = np.arange(np.size(frequencies))[mess_formfactors > 0]
    kk_phase = np.zeros_like(frequencies)
    for i in np.arange(np.size(indices)):
        #avouid 0 in nenner
        ff_here = np.delete(these_ff, i)
        freqs_here = np.delete(these_freqs, i)
        to_integrate = (np.log(ff_here) - np.log(these_ff[i])) / (these_freqs[i] ** 2 - freqs_here ** 2)
        kk_phase[indices[i]] = 2 * these_freqs[i] / np.pi * np.trapz(to_integrate, x=freqs_here)

    return kk_phase


def model_start(frequencies, model_ff, mess_formfactors, mess_formfactors_noise):
    """
    Starts the Gerchberg Saxton algorithm with a complex ff and gives the final time profil
    """

    zeit = np.fft.fftfreq(2 * np.size(frequencies) - 1, d=frequencies[1] - frequencies[0])
    zeit = np.sort(zeit)
    #tmax = 1/(frequencies[0]-frequencies[1])
    #zeit = tmax*np.linspace(-1/2,1/2,np.size(frequencies))
    start_phase = np.angle(model_ff)
    time_domain = gerchberg_saxton(mess_formfactors, mess_formfactors_noise, start_phase)
    return zeit, time_domain


def norm_current(zeit, profil, charge):
    """
    Normalizes the longitudinal profile to the actual charge
    """
    return profil / np.trapz(profil, x=zeit) * charge


def cleanup_formfactor(frequencies, formfactors, formfactors_noise, formfactors_det_limit,
                       channels_to_remove=[104, 135], wanted_time_res=2e-15, wanted_time_frame=2e-12,
                       high_inter_last=1, model_last_index=-1, smooth_window=9):
    """
    These function cleans up the form factor and interpolates it in frequency
    """
    frequencies, formfactors, formfactors_noise, formfactors_det_limit = remove_certain_channels(channels_to_remove,
                                                                                                 [frequencies,
                                                                                                  formfactors,
                                                                                                  formfactors_noise,
                                                                                                  formfactors_det_limit])
    sort_freq, sorted_stuff = sort_data(frequencies, [formfactors, formfactors_noise, formfactors_det_limit])
    sign_frequencies, sign_formfactors, sign_formfactors_noise = remove_loners(sort_freq, sorted_stuff[0],
                                                                               sorted_stuff[1], sorted_stuff[2])
    #New FrequencyGrid
    freq_d = 1 / wanted_time_frame
    freq_max = 1 / wanted_time_res
    n_freqs = int(freq_max / freq_d)
    #Frequencies must contain 0!
    new_freqs = np.linspace(0, freq_max, n_freqs)
    #Extrapolations
    extra_low_freqs, extra_low_ff, extra_low_ff_noise = extrapolation_low(new_freqs, sign_frequencies, sign_formfactors,
                                                                          sign_formfactors_noise)
    extra_high_freqs, extra_high_ff, extra_high_ff_noise = extrapolation_high(new_freqs, extra_low_freqs, extra_low_ff,
                                                                              extra_low_ff_noise,
                                                                              index_to_cut=high_inter_last)
    #Smoothing & Interpolation
    smooth_ff = smooth_formfactors(extra_high_ff, window_size=smooth_window)
    final_ff, final_ff_noise = make_linear_spacing(new_freqs, extra_high_freqs, smooth_ff, extra_high_ff_noise)
    return new_freqs, final_ff, final_ff_noise


def master_recon(frequencies, formfactors, formfactors_noise, formfactors_det_limit, charge, method='KKstart',
                 channels_to_remove=[104, 135], wanted_time_res=2e-15, wanted_time_frame=2e-12,
                 high_inter_last=1, model_last_index=-1, smooth_window=9, phase_noise_start=None, show_plots=True):
    """
    Master Function to reconstruct the current profile from a given formfactor
    measurement. Returns time and current profile
    """

    frequencies, formfactors, formfactors_noise, formfactors_det_limit = remove_certain_channels(channels_to_remove,
                                                                                                 [frequencies,
                                                                                                  formfactors,
                                                                                                  formfactors_noise,
                                                                                                  formfactors_det_limit])
    sort_freq, sorted_stuff = sort_data(frequencies, [formfactors, formfactors_noise, formfactors_det_limit])
    sign_frequencies, sign_formfactors, sign_formfactors_noise = remove_loners(sort_freq, sorted_stuff[0],
                                                                               sorted_stuff[1], sorted_stuff[2])
    #New FrequencyGrid
    freq_d = 1 / wanted_time_frame
    freq_max = 1 / wanted_time_res
    n_freqs = int(freq_max / freq_d)
    #Frequencies must contain 0!
    new_freqs = np.linspace(0, freq_max, n_freqs)
    #Extrapolations
    extra_low_freqs, extra_low_ff, extra_low_ff_noise = extrapolation_low(new_freqs, sign_frequencies, sign_formfactors,
                                                                          sign_formfactors_noise)
    extra_high_freqs, extra_high_ff, extra_high_ff_noise = extrapolation_high(new_freqs, extra_low_freqs, extra_low_ff,
                                                                              extra_low_ff_noise,
                                                                              index_to_cut=high_inter_last)
    #Smoothing & Interpolation
    smooth_ff = smooth_formfactors(extra_high_ff, window_size=smooth_window)
    final_ff, final_ff_noise = make_linear_spacing(new_freqs, extra_high_freqs, smooth_ff, extra_high_ff_noise)
    if method == 'modell':
        #Modelfit
        start_ff = find_best_modelff(modellist, new_freqs, final_ff, sign_frequencies[model_last_index])
        #return start_ff
    elif method == 'KKstart':
        #print('Drin')
        kk_phase = kramers_kronig_phase(new_freqs, final_ff)
        start_ff = final_ff * np.exp(kk_phase * complex(0, 1))
        #return start_ff
    else:
        print('Method not defined')
    if phase_noise_start != None:
        #alle start phasen ueber 20THz random verteilen
        crit_freq = phase_noise_start
        n_phases = np.size(new_freqs[new_freqs > crit_freq])
        rand_phases = np.random.rand(n_phases) * 2 * np.pi
        start_ff[new_freqs > crit_freq] = np.abs(start_ff[new_freqs > crit_freq]) * np.exp(complex(0, 1) * rand_phases)

    #Reconstruction
    recon_time, recon_prof = model_start(new_freqs, start_ff, final_ff, final_ff_noise)
    #Normalization to current
    current = norm_current(recon_time, recon_prof, charge)

    #get formfactor of reconstruction
    #wx, recon_ff = nils_fft(recon_time,current/charge)
    wx, recon_ff = nils_fft(np.arange(np.size(recon_time)), recon_prof)
    recon_ff = recon_ff * np.sqrt(2 * np.pi)
    recon_ff = recon_ff[:int((np.size(recon_ff) + 1) / 2)]
    recon_ff = np.abs(recon_ff)

    #RMS Time Value
    t_rms = np.sqrt(np.abs(second_moment(current, recon_time)))
    if show_plots == True:
        fig, axes = plt.subplots(nrows=2, ncols=1)
        #Raw Data
        axes[0].errorbar(sort_freq * 1e-12, sorted_stuff[0], yerr=sorted_stuff[1], fmt='o', label='Raw Data', color='b')
        axes[0].fill_between(sort_freq * 1e-12, 1e-4, sorted_stuff[2], label='Det Limit', color='gray', alpha=0.3)
        #Before Recon
        axes[0].errorbar(new_freqs * 1e-12, final_ff, yerr=final_ff_noise, fmt='-', label='Input to Recon.',
                         color='orange')
        #Modell
        if method == 'modell':
            axes[0].plot(new_freqs * 1e-12, np.abs(start_ff), label='Model', color='g')
        #Results of Recon
        axes[0].plot(new_freqs * 1e-12, recon_ff, '--', label='Result of Recon.', color='r')
        #To look nice
        axes[0].set_xscale('log')
        axes[0].set_yscale('log')
        axes[0].set_xlim([0.5, 100])
        axes[0].set_xlabel('f (THz)')
        axes[0].set_ylabel('|F|')
        axes[0].set_ylim([5e-3, 2])
        axes[0].legend(loc='lower left')

        #Time Domain
        axes[1].plot(recon_time * 1e15, current * 1e-3, label=r'$\sigma_t = $' + str(int(round(t_rms * 1e15))) + ' fs')
        axes[1].set_xlabel('t (fs)')
        axes[1].set_ylabel('I (kA)')
        axes[1].set_xlim([-10 * t_rms * 1e15, 10 * t_rms * 1e15])
        #axes[1].set_xlim([-200,200])
        axes[1].legend()
        fig.tight_layout()

    return recon_time, current, t_rms


"""
simulate spectrometer

Function to simulate the measured formfactor of the spectrometer for input of
s (m) and I (amps)
"""


def get_crisp_signal(s, current, n_shots=1, which_set='high'):
    """Simulates the signal of the spectrometer of a single shot including the
    noise of the spectrometer.

    Parameters
    ----------
    s :         numpy array of floats
                longitudinal bunch coordinate.
    current :   numpy array of floats with same length as s
                current in ampere
    n_shots:    int, optional
                number of shots to average from. Reduces spectrometer noise.
    which_set:  str, optional
                Specifies the grating set of the spectrometer to use. Keep in
                mind that experimentally you can only get one at a time.
                Keywords are: 'high', 'low' and 'both'

    Returns
    ------
    results :   float (array)
                resulting electric field at distance r
    """
    # center s so that 0 is in the middle
    s_centered = s - np.sum(s * current) / np.sum(current)
    # get charge
    t_in = s_centered / 3e8
    charge = np.trapz(current, x=t_in)
    # linear time spacing
    t_spec = np.arange(-2, 2, 0.001) * 1e-12
    # needs to be sorted
    current = current[np.argsort(t_in)]
    t_in = np.sort(t_in)
    cur_spec = np.interp(t_spec, t_in, current, left=0, right=0)
    # normalize to one
    cur_spec = cur_spec / np.trapz(cur_spec, x=t_spec)

    # plt.figure('Normalized and intrapolated')
    # plt.plot(t_spec*1e15,cur_spec*1e-3)

    # fft
    comp_ff = np.fft.fft(cur_spec)
    scale_factor = (t_spec[1] - t_spec[0])
    comp_ff = comp_ff * scale_factor
    freqs = np.fft.fftfreq(np.size(cur_spec), t_spec[1] - t_spec[0])
    comp_ff = comp_ff[freqs >= 0]
    freqs = freqs[freqs >= 0]

    # plt.figure()
    # plt.plot(freqs*1e-12,np.abs(comp_ff))
    # plt.xscale('log')
    # plt.yscale('log')

    # interpolate to spectrometer range
    freqs_spectrometer = np.array([6.84283010e+11, 6.86342276e+11, 6.89291238e+11, 6.93231268e+11,
                                   6.97867818e+11, 7.03183896e+11, 7.09148197e+11, 7.15734438e+11,
                                   7.23107481e+11, 7.31240520e+11, 7.40197234e+11, 7.50015830e+11,
                                   7.60740352e+11, 7.72411055e+11, 7.85016183e+11, 7.98657943e+11,
                                   8.13405379e+11, 8.29377459e+11, 8.46469777e+11, 8.65200287e+11,
                                   8.85430610e+11, 9.07299763e+11, 9.31292658e+11, 9.57343678e+11,
                                   9.85437362e+11, 1.01603519e+12, 1.04944212e+12, 1.08594355e+12,
                                   1.12511418e+12, 1.16825162e+12, 1.25452235e+12, 1.25941010e+12,
                                   1.26511563e+12, 1.27198544e+12, 1.28040955e+12, 1.29006363e+12,
                                   1.30097204e+12, 1.31295400e+12, 1.32651050e+12, 1.34139271e+12,
                                   1.35760346e+12, 1.37547358e+12, 1.39503008e+12, 1.41656334e+12,
                                   1.43974321e+12, 1.46457178e+12, 1.49167731e+12, 1.52101420e+12,
                                   1.55222847e+12, 1.58655331e+12, 1.62359096e+12, 1.66369362e+12,
                                   1.70718601e+12, 1.75508152e+12, 1.80652971e+12, 1.86252358e+12,
                                   1.92362000e+12, 1.99001009e+12, 2.06242973e+12, 2.14196419e+12,
                                   2.28119024e+12, 2.28903842e+12, 2.29950939e+12, 2.31261258e+12,
                                   2.32772918e+12, 2.34515972e+12, 2.36498904e+12, 2.38681143e+12,
                                   2.41113054e+12, 2.43827617e+12, 2.46831212e+12, 2.50106736e+12,
                                   2.53639935e+12, 2.57499338e+12, 2.61694409e+12, 2.66265777e+12,
                                   2.71189970e+12, 2.76495993e+12, 2.82181770e+12, 2.88430880e+12,
                                   2.95157763e+12, 3.02570939e+12, 3.10515027e+12, 3.19080363e+12,
                                   3.28411150e+12, 3.38505628e+12, 3.49557917e+12, 3.61569473e+12,
                                   3.74494818e+12, 3.87221178e+12, 3.87605079e+12, 3.89182771e+12,
                                   3.91112792e+12, 3.93302151e+12, 3.95880137e+12, 3.98757241e+12,
                                   4.02090481e+12, 4.05868059e+12, 4.09887219e+12, 4.14427297e+12,
                                   4.19590008e+12, 4.25159913e+12, 4.31184945e+12, 4.37798056e+12,
                                   4.45017789e+12, 4.52767991e+12, 4.61054903e+12, 4.70100308e+12,
                                   4.79815306e+12, 4.90363011e+12, 5.01847576e+12, 5.14289228e+12,
                                   5.27853835e+12, 5.42497233e+12, 5.58394521e+12, 5.75566213e+12,
                                   5.94270061e+12, 6.14675659e+12, 6.36808841e+12, 6.60398418e+12,
                                   6.84239080e+12, 6.86645607e+12, 6.89844358e+12, 6.93795435e+12,
                                   6.98381615e+12, 7.03665301e+12, 7.09581595e+12, 7.16182658e+12,
                                   7.23558808e+12, 7.31656935e+12, 7.40588391e+12, 7.50366620e+12,
                                   7.61017924e+12, 7.72725265e+12, 7.85313450e+12, 7.98952729e+12,
                                   8.13636493e+12, 8.29743559e+12, 8.46741872e+12, 8.65427554e+12,
                                   8.85635877e+12, 9.07462058e+12, 9.31463114e+12, 9.57583282e+12,
                                   9.85392517e+12, 1.01581001e+13, 1.04872716e+13, 1.08470040e+13,
                                   1.12261054e+13, 1.13980414e+13, 1.14448608e+13, 1.15007118e+13,
                                   1.15668861e+13, 1.16436805e+13, 1.16745582e+13, 1.17301474e+13,
                                   1.18286343e+13, 1.19375252e+13, 1.20584689e+13, 1.21926665e+13,
                                   1.23420542e+13, 1.25037214e+13, 1.26827481e+13, 1.28758007e+13,
                                   1.30836435e+13, 1.33098200e+13, 1.35576821e+13, 1.38245741e+13,
                                   1.41103234e+13, 1.44217346e+13, 1.47624195e+13, 1.51316322e+13,
                                   1.55258557e+13, 1.59552050e+13, 1.64196917e+13, 1.69285923e+13,
                                   1.74782780e+13, 1.80802622e+13, 1.87316513e+13, 1.94415752e+13,
                                   2.05304449e+13, 2.06041596e+13, 2.06980683e+13, 2.08140040e+13,
                                   2.09510109e+13, 2.11093505e+13, 2.12864874e+13, 2.14804199e+13,
                                   2.17008944e+13, 2.19489705e+13, 2.22216055e+13, 2.25161992e+13,
                                   2.28301018e+13, 2.31801788e+13, 2.35588868e+13, 2.39691136e+13,
                                   2.44111492e+13, 2.48897567e+13, 2.54005500e+13, 2.59618741e+13,
                                   2.65693466e+13, 2.72227653e+13, 2.79406515e+13, 2.87222610e+13,
                                   2.95520370e+13, 3.04799241e+13, 3.14692524e+13, 3.25463980e+13,
                                   3.36766529e+13, 3.41763020e+13, 3.43153008e+13, 3.44839431e+13,
                                   3.46803032e+13, 3.49137251e+13, 3.49825754e+13, 3.51760204e+13,
                                   3.54682908e+13, 3.57958567e+13, 3.61610300e+13, 3.65661507e+13,
                                   3.70094847e+13, 3.74942425e+13, 3.80275217e+13, 3.86055890e+13,
                                   3.92357545e+13, 3.99218111e+13, 4.06583840e+13, 4.14533934e+13,
                                   4.23052569e+13, 4.32362878e+13, 4.42368558e+13, 4.53735986e+13,
                                   4.65444168e+13, 4.78371060e+13, 4.92264177e+13, 5.07656062e+13,
                                   5.24232320e+13, 5.41042876e+13, 5.61221152e+13, 5.82673400e+13])

    ff_spectrometer = np.interp(freqs_spectrometer, freqs, np.abs(comp_ff))

    # plt.figure('Spectrometer Signal')
    # plt.plot(freqs_spectrometer*1e-12, ff_spectrometer, 'o')
    # plt.xscale('log')
    # plt.yscale('log')

    # add spectrometer noise
    # response in V/nC²
    spec_response = np.array([2.27139720e-02, 2.63574891e-02, 1.91519673e-02, 2.68880698e-02,
                              2.96823444e-02, 2.50777417e-02, 5.14215877e-02, 3.65048653e-02,
                              5.01399530e-02, 4.92891672e-02, 6.28404705e-02, 5.90190272e-02,
                              7.09957402e-02, 9.88634708e-02, 1.14630994e-01, 1.23265766e-02,
                              1.23461331e-01, 1.36320335e-01, 1.25770119e-01, 1.32215392e-01,
                              1.36440279e-01, 1.34883965e-01, 1.31848572e-01, 1.64531877e-01,
                              1.55915488e-01, 1.92955737e-01, 2.10324614e-01, 2.06284403e-01,
                              1.21629041e-01, 1.28159124e-01, 2.82388474e-02, 3.24896486e-02,
                              6.32553505e-02, 7.70183540e-02, 8.65472943e-02, 1.15282182e-01,
                              1.20753130e-01, 1.58432901e-01, 1.75611721e-01, 1.91465810e-01,
                              2.57557223e-01, 2.96939075e-01, 3.14254548e-01, 3.66881036e-01,
                              3.91964519e-01, 3.94565181e-01, 3.35777982e-01, 4.34468414e-01,
                              4.34228563e-01, 4.55406870e-01, 6.28887904e-01, 6.28342309e-01,
                              7.46633419e-01, 7.66892682e-01, 7.07245311e-01, 6.48177114e-01,
                              4.44611893e-01, 5.51623941e-01, 7.26412461e-01, 6.13676906e-01,
                              6.79542586e-02, 8.16807490e-02, 1.32899393e-01, 1.91507017e-01,
                              1.90696482e-01, 2.00147063e-01, 2.73266787e-01, 2.51203348e-01,
                              2.91641304e-01, 3.02340215e-01, 2.90033553e-01, 3.56812708e-01,
                              3.80332636e-01, 4.71988270e-01, 5.71053105e-01, 5.76061366e-01,
                              7.42794124e-01, 8.28076940e-01, 8.61271616e-01, 9.34415109e-01,
                              9.88661829e-01, 1.29302949e+00, 1.68370088e+00, 1.60985763e+00,
                              1.66388597e+00, 2.18810117e+00, 1.77063413e+00, 2.22234386e+00,
                              1.79544513e+00, 6.49211052e-01, 1.12905733e-01, 1.35219749e-01,
                              1.77552881e-01, 2.33864266e-01, 2.76136407e-01, 3.08371345e-01,
                              4.22577865e-01, 4.50832291e-01, 4.07343365e-01, 4.54518985e-01,
                              4.74195250e-01, 5.14916496e-01, 5.82809637e-01, 6.85511764e-01,
                              7.85862983e-01, 9.18192559e-01, 1.08807514e+00, 1.32820958e+00,
                              1.63915706e+00, 1.88506672e+00, 2.11140637e+00, 2.56819078e+00,
                              3.47029067e+00, 4.43909591e+00, 4.59933654e+00, 5.20503438e+00,
                              4.63311775e+00, 4.98549453e+00, 6.05139362e+00, 2.00813681e+00,
                              3.17399571e-01, 4.80405496e-01, 6.33070146e-01, 7.01504333e-01,
                              8.78360578e-01, 1.04252778e+00, 1.33420246e+00, 1.44414774e+00,
                              1.60301216e+00, 1.95838326e+00, 1.97299122e+00, 2.01505983e+00,
                              2.48201305e+00, 3.35207178e+00, 4.50221954e+00, 1.02154423e-01,
                              5.83456121e+00, 6.21937385e+00, 6.72302208e+00, 6.93572009e+00,
                              7.72672143e+00, 9.53759783e+00, 9.70607794e+00, 9.88128478e+00,
                              1.02574353e+01, 1.17366859e+01, 9.75250128e+00, 8.49536317e+00,
                              5.87634647e+00, 6.90542536e-01, 8.84953809e-01, 1.10143838e+00,
                              1.32548040e+00, 1.39920693e+00, 2.64741690e-01, 1.49174256e+00,
                              1.42851071e+00, 1.61623263e+00, 2.14723132e+00, 2.49236663e+00,
                              2.41236131e+00, 2.52993530e+00, 1.95052864e+00, 2.08077797e+00,
                              2.46425871e+00, 1.86427999e+00, 1.29883560e+00, 1.61806660e+00,
                              2.69349944e+00, 5.26651344e+00, 8.60076710e+00, 1.23955411e+01,
                              1.55537653e+01, 1.57520786e+01, 1.55318570e+01, 1.63010277e+01,
                              1.95997146e+01, 1.89657798e+01, 1.98173958e+01, 1.28444056e+01,
                              1.05663513e+00, 1.73837417e+00, 2.38571066e+00, 3.12509262e+00,
                              5.16370343e+00, 6.27279615e+00, 4.55271819e+00, 2.30736012e+00,
                              1.00261293e+00, 1.28871366e+00, 6.82622230e+00, 1.52316196e+01,
                              1.30645368e+01, 1.56890352e+01, 1.39606414e+01, 1.52718010e+01,
                              7.30423998e+00, 1.86593830e+01, 2.01944545e+01, 2.18768013e+01,
                              2.09667009e+01, 2.15890054e+01, 2.39164749e+01, 2.48425545e+01,
                              2.30255779e+01, 1.71761298e+01, 2.33330662e+01, 2.32725960e+01,
                              2.10340416e+01, 5.24993968e-01, 1.63131206e+00, 3.44781069e+00,
                              2.29974095e+00, 4.99462286e+00, 2.32302429e+00, 4.96744300e+00,
                              2.68405786e+00, 8.62683761e+00, 1.15968873e+01, 1.43637056e+01,
                              1.16110817e+01, 9.31426757e+00, 1.07818775e+01, 1.02543508e+01,
                              7.20074507e+00, 1.01373496e+01, 1.66319551e+01, 5.87900394e+00,
                              1.78322017e+01, 6.29817365e+00, 9.87433101e-01, 8.83113064e+00,
                              6.72258056e+00, 7.85860449e+00, 9.53416250e+00, 2.59286607e+01,
                              2.66550560e+01, 2.05970702e+01, 2.45164563e+00, 1.38458228e+00])

    # formfactor to adc-signal
    adc_sig = ff_spectrometer ** 2 * (charge * 1e9) ** 2 * spec_response  # in V
    # plt.figure('ADC Signal')
    # plt.plot(adc_sig)
    # add_noise
    elec_noise = 1.2e-3  # V
    elec_noise = elec_noise / np.sqrt(n_shots)
    adc_noise = np.random.randn(np.size(ff_spectrometer)) * elec_noise
    adc_total = adc_sig + adc_noise

    # And back to formfactor
    final_ff = 1 / (charge * 1e9) * np.sqrt(np.abs(adc_total) / spec_response) * np.sign(adc_total)
    # noise on form factor
    ff_noise = 0.5 * final_ff * elec_noise * np.sqrt(n_shots) / adc_total  # fehlerfortpflanzung
    det_lim = 1 / (charge * 1e9) * np.sqrt(np.abs(elec_noise) / spec_response)  # noise floor
    # plt.figure()
    # plt.plot(freqs_spectrometer*1e-12, spec_response, 'o')
    # plt.xscale('log')
    # plt.yscale('log')

    # plt.figure('Spectrometer Signal')
    # plt.plot(freqs_spectrometer*1e-12, final_ff, 'o')
    # plt.xscale('log')
    # plt.yscale('log')
    if which_set == 'low':
        freqs_to_return = freqs_spectrometer[:120]
        ff_to_return = final_ff[:120]
        det_lim = det_lim[:120]
        ff_noise = ff_noise[:120]
    elif which_set == 'high':
        freqs_to_return = freqs_spectrometer[120:]
        ff_to_return = final_ff[120:]
        det_lim = det_lim[120:]
        ff_noise = ff_noise[120:]
    elif which_set == 'both':
        freqs_to_return = freqs_spectrometer
        ff_to_return = final_ff
        det_lim = det_lim
        ff_noise = ff_noise
    else:
        print('Keyword not known!')

    return np.array([freqs_to_return, ff_to_return, ff_noise, det_lim])


