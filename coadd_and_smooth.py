from __future__ import division

from astropy.io import fits
import numpy as np
from astropy.convolution import convolve, Gaussian1DKernel

import os
import sys

import matplotlib.pyplot as plt
import plot_final_spectra as pfs

def get_spectra(obj_hdu):

    jsplt = obj_hdu[0].header["WAT2_001"].split('=')[-1].split(' ')
    hsplt = obj_hdu[0].header["WAT2_002"].split('=')[-1].split(' ')
    ksplt = obj_hdu[0].header["WAT2_003"].split('=')[-1].split(' ')

    jstart = float(jsplt[4])
    delta_j = float(jsplt[5])

    hstart = float(hsplt[4])
    delta_h = float(hsplt[5])

    kstart = float(ksplt[4])
    delta_k = float(ksplt[5])

    num = 2048

    jwav = [jstart+i*delta_j for i in range(2048)]
    jspec = obj_hdu[0].data[0,0]

    hwav = [hstart+i*delta_h for i in range(2048)]
    hspec = obj_hdu[0].data[0,1]

    kwav = [kstart+i*delta_k for i in range(2048)]
    kspec = obj_hdu[0].data[0,2]

    return jspec, hspec, kspec, np.asarray(jwav), np.asarray(hwav), np.asarray(kwav), jstart, hstart, kstart

def smoothspectra(spec, width=2):

    gauss_kernel = Gaussian1DKernel(width)  # width=1 is the same as no smoothing
    smoothed_spec = convolve(spec, gauss_kernel)

    return smoothed_spec

def coadd(spec_ab, spec_ba, wav_ab, wav_ba, wav_grid, avg_delta):

    coadded_spec = np.zeros(len(wav_grid))

    for u in range(len(wav_grid)):

        idx_ab = np.where((wav_ab >= wav_grid[u] - avg_delta/2) & (wav_ab < wav_grid[u] + avg_delta/2))[0]
        idx_ba = np.where((wav_ba >= wav_grid[u] - avg_delta/2) & (wav_ba < wav_grid[u] + avg_delta/2))[0]

        conc_arr = np.concatenate((spec_ab[idx_ab], spec_ba[idx_ba]), axis=0)
        coadded_spec[u] = np.median(conc_arr)

    return coadded_spec

if __name__ == '__main__':

    # read in fits file and get wav calibrated spectra
    # give it the filename which has the dispersion corrected spectra
    obj_name = 'xl53'
    slitpos = 'AB'
    obj_filename = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/20170511/' + obj_name + '_' \
    + slitpos + '_tellinterp_dispcor.fits'
    obj_spec_ab = fits.open(obj_filename)

    slitpos = 'BA'
    obj_filename = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/20170511/' + obj_name + '_' \
    + slitpos + '_tellinterp_dispcor.fits'
    obj_spec_ba = fits.open(obj_filename)

    jspec_ab, hspec_ab, kspec_ab, jwav_ab, hwav_ab, kwav_ab, jstart_ab, hstart_ab, kstart_ab = get_spectra(obj_spec_ab)
    jspec_ba, hspec_ba, kspec_ba, jwav_ba, hwav_ba, kwav_ba, jstart_ba, hstart_ba, kstart_ba = get_spectra(obj_spec_ba)

    # Now coadd the A and B position spectra
    # smooth the spectra and then coadd (I'm not sure if this is the correct order of doing things.)
    jspec_ab_smooth = smoothspectra(jspec_ab)
    hspec_ab_smooth = smoothspectra(hspec_ab)
    kspec_ab_smooth = smoothspectra(kspec_ab)

    jspec_ba_smooth = smoothspectra(jspec_ba)
    hspec_ba_smooth = smoothspectra(hspec_ba)
    kspec_ba_smooth = smoothspectra(kspec_ba)

    # create a lambda arrays that are consistent for both
    avg_jstart = int( min(jstart_ab, jstart_ba) - 1)
    avg_hstart = int( min(hstart_ab, hstart_ba) - 1)
    avg_kstart = int( min(kstart_ab, kstart_ba) - 1)

    avg_delta = 6.0

    num_j = int( (14837 - avg_jstart) / avg_delta )
    num_h = int( (18520 - avg_hstart) / avg_delta )
    num_k = int( (24658 - avg_kstart) / avg_delta )

    jwav_grid = [avg_jstart + i*avg_delta for i in range(num_j)]
    hwav_grid = [avg_hstart + i*avg_delta for i in range(num_h)]
    kwav_grid = [avg_kstart + i*avg_delta for i in range(num_k)]

    # now move along the wavelength arrays and coadd spectra
    j_coadd = coadd(jspec_ab, jspec_ba, jwav_ab, jwav_ba, jwav_grid, avg_delta)
    h_coadd = coadd(hspec_ab, hspec_ba, hwav_ab, hwav_ba, hwav_grid, avg_delta)
    k_coadd = coadd(kspec_ab, kspec_ba, kwav_ab, kwav_ba, kwav_grid, avg_delta)

    pfs.plotspec(jwav_grid, j_coadd, obj_name, 'j', 'coadd')
    pfs.plotspec(hwav_grid, h_coadd, obj_name, 'h', 'coadd')
    pfs.plotspec(kwav_grid, k_coadd, obj_name, 'k', 'coadd')

    sys.exit(0)