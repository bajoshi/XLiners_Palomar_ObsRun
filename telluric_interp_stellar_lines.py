from __future__ import division

from astropy.io import fits
import numpy as np
from scipy.interpolate import interp1d

import os
import sys
import copy

import matplotlib.pyplot as plt 

import rm_nan_norm as rnn

def interp_and_plot(line_left, line_right, tell_spec_for_line):

    line_x = np.concatenate((line_left, line_right), axis=0)
    line_y = tell_spec_for_line[line_x]
    line_func = interp1d(line_x, line_y, kind='linear')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    x_intp = np.arange(line_left[0], line_right[-1])
    ax.plot(x_intp, tell_spec_for_line[line_left[0]: line_right[-1]])
    ax.plot(x_intp, line_func(x_intp))

    plt.show()

    return line_func, x_intp

def plot_tell_after_interp(tell_old, tell_new):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    xvals = np.arange(len(tell_old))
    ax.plot(xvals, tell_old, 'b')
    ax.plot(xvals, tell_new, 'g')

    plt.show()

    return None

def interpolate_over_stellar_lines(tell_hdu):
    """
    This function will interpolate over H recomb absorption lines seen in the spectrum of 
    the telluric standard star. 

    We only want to divide out the telluric absorption features, 
    which come from the earth's atmosphere, in the galaxy spectrum using the telluric 
    standard star. But these absorption features in the standard star come from the stellar 
    atmosphere, NOT the earth's atmosphere. Therefore, we are interpolating over them.
    i.e. our assumption of a flat spectrum for the standard star was incorrect; which we
    knew but we didn't know that the absorption features were deep enough to cause 
    artificial emission lines to appear in the final galaxy spectrum after the telluric
    division.

    It expects you to give it the J, H, and K band spectra for the telluric standard star.
    """

    # get the normalized spectra first and copy (see python docs on deep copy) 
    # them for plotting purposes
    tell_spec_j, tell_spec_h, tell_spec_k = rnn.insert_nan_norm_tell(tell_hdu)

    tell_spec_j_old = copy.deepcopy(tell_spec_j)
    tell_spec_h_old = copy.deepcopy(tell_spec_h)
    tell_spec_k_old = copy.deepcopy(tell_spec_k)

    # identify stellar lines (and corresponding pixels) that you want to interpolate over.
    # since we are working here with spectra that are not wavelength 
    # calibrated, the slicing will be done using pixels.
    # I think the slicing pixels should be the same for all stars but
    # this should be confirmed by eye for each star that this code is used for.
    # Info from:
    # http://www.not.iac.es/instruments/notcam/ReferenceInfo/  # home page
    # http://www.not.iac.es/instruments/notcam/ReferenceInfo/recomb_menu.html  # atomic recombination lines
    # http://www.not.iac.es/instruments/notcam/ReferenceInfo/h2_lines.html  # molecular hydrogen lines
    # although we certainly don't see any molecular Hydrogen lines from the white dwarfs 
    # atmosphere (not surprising... lol), I'm keeping this page here because I need the wavelengths
    # for hte H2 lines when identifying lines in the galaxy spectra.
    # Variables names:
    # Br -- Brackett line
    # Pa -- Paschen line
    # left and right indicate numpy array indices that should be used for interpolation
    # on either side of the absorption feature.

    #Pa_alpha =  # rest wav = 1.8756 microns

    Pa_beta_left = np.arange(1120,1140) # rest wav = 1.2822 microns  # J band
    Pa_beta_right = np.arange(1190,1220)

    Br_gamma_left = np.arange(990,1021) # rest wav = 2.1661 microns  # K band
    Br_gamma_right = np.arange(1070,1080)

    Br_zeta_left = np.arange(505,509) # rest wav = 1.7367 microns  # H band
    Br_zeta_right = np.arange(551,553)

    Br_eta_left = np.arange(729,759) # rest wav = 1.6811 microns  # H band
    Br_eta_right = np.arange(823,836)

    Br_theta_left = np.arange(951,955) # rest wav = 1.6412 microns  # H band
    Br_theta_right = np.arange(1010,1060)

    Br_iota_left = np.arange(1020,1060) # rest wav = 1.6114 microns  # H band
    Br_iota_right = np.arange(1144,1145)

    Br_kappa_left = np.arange(1186,1196) # rest wav = 1.5885 microns  # H band
    Br_kappa_right = np.arange(1245,1256)

    Br_lambda_left = np.arange(1247,1257) # rest wav = 1.5705 microns  # H band
    Br_lambda_right = np.arange(1331,1347)

    Br_mu_left = np.arange(1331,1350) # rest wav = 1.5561 microns  # H band
    Br_mu_right = np.arange(1399,1415)

    # interpolate and plot the result of interpolation
    Pa_beta_linefunc, Pa_beta_xvals = interp_and_plot(Pa_beta_left, Pa_beta_right, tell_spec_j)

    Br_gamma_linefunc, Br_gamma_xvals = interp_and_plot(Br_gamma_left, Br_gamma_right, tell_spec_k)

    Br_zeta_linefunc, Br_zeta_xvals = interp_and_plot(Br_zeta_left, Br_zeta_right, tell_spec_h)
    Br_eta_linefunc, Br_eta_xvals = interp_and_plot(Br_eta_left, Br_eta_right, tell_spec_h)
    Br_theta_linefunc, Br_theta_xvals = interp_and_plot(Br_theta_left, Br_theta_right, tell_spec_h)
    Br_iota_linefunc, Br_iota_xvals = interp_and_plot(Br_iota_left, Br_iota_right, tell_spec_h)
    Br_kappa_linefunc, Br_kappa_xvals = interp_and_plot(Br_kappa_left, Br_kappa_right, tell_spec_h)
    Br_lambda_linefunc, Br_lambda_xvals = interp_and_plot(Br_lambda_left, Br_lambda_right, tell_spec_h)
    Br_mu_linefunc, Br_mu_xvals = interp_and_plot(Br_mu_left, Br_mu_right, tell_spec_h)

    # stitch back the interpolated section
    tell_spec_j[Pa_beta_xvals] = Pa_beta_linefunc(Pa_beta_xvals)

    tell_spec_k[Br_gamma_xvals] = Br_gamma_linefunc(Br_gamma_xvals)

    tell_spec_h[Br_zeta_xvals] = Br_zeta_linefunc(Br_zeta_xvals)
    tell_spec_h[Br_eta_xvals] = Br_eta_linefunc(Br_eta_xvals)
    tell_spec_h[Br_theta_xvals] = Br_theta_linefunc(Br_theta_xvals)
    tell_spec_h[Br_iota_xvals] = Br_iota_linefunc(Br_iota_xvals)
    tell_spec_h[Br_kappa_xvals] = Br_kappa_linefunc(Br_kappa_xvals)
    tell_spec_h[Br_lambda_xvals] = Br_lambda_linefunc(Br_lambda_xvals)
    tell_spec_h[Br_mu_xvals] = Br_mu_linefunc(Br_mu_xvals)

    plot_tell_after_interp(tell_spec_j_old, tell_spec_j)
    plot_tell_after_interp(tell_spec_h_old, tell_spec_h)
    plot_tell_after_interp(tell_spec_k_old, tell_spec_k)

    tell_hdu[0].data[0,0] = tell_spec_j
    tell_hdu[0].data[0,1] = tell_spec_h
    tell_hdu[0].data[0,2] = tell_spec_k

    return tell_hdu

def do_all(work_dir, tell_name, slitpos):

    tell_filename = work_dir + tell_name + '_' + slitpos + '_bksub.fits'
    tell_star_spec = fits.open(tell_filename)

    tell_star_spec = interpolate_over_stellar_lines(tell_star_spec)
    tell_star_spec.writeto(tell_filename.replace('_bksub.fits', '_stellarlines_interp_norm_bksub.fits'), clobber=True)
    tell_star_spec.close()
    del tell_star_spec

    return None

if __name__ == '__main__':

    slitpos = 'BA'
    ext_dir = '/Volumes/Bhavins_backup/ipac/Palomar_data/test/'
    date = ''
    
    # ------------- first telluric --------------- #
    tell_name = 'hip64248'
    do_all(ext_dir, date, tell_name, slitpos)

    # ------------- second telluric --------------- #
    #tell_name = 'hip65280'
    #do_all(tell_name, slitpos)

    sys.exit(0)