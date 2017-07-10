from __future__ import division

from astropy.io import fits
import numpy as np

import os
import sys

import matplotlib.pyplot as plt 

def plotspec(spec):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.arange(len(spec)), spec)
    plt.show()

    return None

def overplot_telluric(tellspec, objspec, poly):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    tellspec *= 10

    xvals = np.arange(len(tellspec))

    ax.plot(xvals, tellspec, color='g')
    ax.plot(xvals, objspec, color='b')

    ax.plot(xvals, objspec/tellspec, color='k', alpha=0.5)
    ax.plot(xvals, objspec * 0.5/poly(xvals), color='r', alpha=0.5)

    #ax.axvline(x=1170, ls='--')

    #ax.set_yscale('log')

    plt.show()

    return None

def plot_tellurics(tellspec1, tellspec2):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(np.arange(len(tellspec1)), tellspec1, color='g') 
    ax.plot(np.arange(len(tellspec2)), tellspec2, color='b') 

    ax.set_yscale('log')

    plt.show()

    return None

def insert_nan_norm_tell(tell_spec_hdu):

    tell_jspec = tell_spec_hdu[0].data[0,0]
    tell_hspec = tell_spec_hdu[0].data[0,1]
    tell_kspec = tell_spec_hdu[0].data[0,2]

    ## remove 0 or negative numbers and replace with NaN
    tell_invalid_jspec = np.where(tell_jspec <= 0)[0]
    tell_invalid_hspec = np.where(tell_hspec <= 0)[0]
    tell_invalid_kspec = np.where(tell_kspec <= 0)[0]

    tell_jspec[tell_invalid_jspec] = np.nan
    tell_hspec[tell_invalid_hspec] = np.nan
    tell_kspec[tell_invalid_kspec] = np.nan

    # normalize star spectrum and write fits file
    tell_jspec_norm = tell_jspec / np.nanmedian(tell_jspec)
    tell_hspec_norm = tell_hspec / np.nanmedian(tell_hspec)
    tell_kspec_norm = tell_kspec / np.nanmedian(tell_kspec)

    tell_spec_hdu[0].data[0,0] = tell_jspec_norm
    tell_spec_hdu[0].data[0,1] = tell_hspec_norm
    tell_spec_hdu[0].data[0,2] = tell_kspec_norm

    return tell_jspec_norm, tell_hspec_norm, tell_kspec_norm

def fit_polynomial_to_spec(spec, order):
    """
    Make sure order is an integer. 
    THis func will return a numpy poly1d object which can simply be used as regular function.

    The presence of NaN's will cause problems with the polynomial fitting,
    so fit the polynomial only to finite values.
    """

    xvals = np.arange(len(spec))

    idx = np.isfinite(spec)

    z = np.polyfit(xvals[idx], spec[idx], order)
    p = np.poly1d(z)

    return p

def plot_poly_fit(spec, poly):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    xvals = np.arange(len(spec))

    ax.plot(xvals, spec, color='g') 
    ax.plot(xvals, poly(xvals), color='b')

    plt.show()

    return None

if __name__ == '__main__':
    
    # ------------- first telluric --------------- #
    tell_name = 'hip64248'
    obj_name = 'xl435'
    slitpos = 'AB'

    tell_filename = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/20170511/' + tell_name + '_' + slitpos + '_bksub.fits'
    obj_filename = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/20170511/' + obj_name + '_' + slitpos + '_bksub.fits'
    obj_spec = fits.open(obj_filename)

    tell_star_spec = fits.open(tell_filename)
    tell_jspec_norm, tell_hspec_norm, tell_kspec_norm = insert_nan_norm_tell(tell_star_spec)

    # ------------- second telluric --------------- #
    tell_name2 = 'hip65280'
    tell_filename2 = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/20170511/' + tell_name2 + '_' + slitpos + '_bksub.fits'
    tell_star_spec2 = fits.open(tell_filename2)

    tell_jspec_norm2, tell_hspec_norm2, tell_kspec_norm2 = insert_nan_norm_tell(tell_star_spec2)

    poly_j = fit_polynomial_to_spec(tell_jspec_norm, 16)
    plot_poly_fit(tell_jspec_norm, poly_j)

    poly_h = fit_polynomial_to_spec(tell_hspec_norm, 15)
    plot_poly_fit(tell_hspec_norm, poly_h)

    poly_k = fit_polynomial_to_spec(tell_kspec_norm, 13)
    plot_poly_fit(tell_kspec_norm, poly_k)

    #tell_star_spec.writeto(tell_filename.replace('.fits', '_norm.fits'), clobber=True)

    # overplot tell and obj to check if the emission lines in the final spectra
    # come from division at wavelengths where telluric absorption exists.
    overplot_telluric(tell_jspec_norm, obj_spec[0].data[0,0], poly_j)
    overplot_telluric(tell_hspec_norm, obj_spec[0].data[0,1], poly_h)
    overplot_telluric(tell_kspec_norm, obj_spec[0].data[0,2], poly_k)

    sys.exit(0)