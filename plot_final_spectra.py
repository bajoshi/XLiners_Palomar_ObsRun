from __future__ import division

from astropy.io import fits
import numpy as np

import os
import sys

import matplotlib.pyplot as plt 

def plotspec(work_dir, wav, spec, obj_name, band, slitpos, redshift):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{\lambda\ [\AA]}$')
    ax.set_ylabel(r'$\mathrm{f\ [arbitrary scale]}$')

    ax.plot(wav, spec, color='k')

    # add vertical lines at the redshifted positions of most likely emission lines
    # define rest wavelengths # in angstroms
    # Only have these for K band right now
    pas_alpha = 18751  # air
    h2_1_0_s3 = 19576
    h2_1_0_s2 = 20338
    h2_1_0_s1 = 21218
    br_gamma = 21655  # air
    br_delta = 19446  # air

    if band == 'j':
        ax.set_xlim(11500, 13500)
        ax.set_ylim(20, 100)
    elif band == 'h':
        ax.set_xlim(14900, 18000)
        ax.set_ylim(100, 500)
    elif band == 'k':
        ax.set_xlim(19000, 24100)
        ax.set_ylim(100, 1000)

        ax.axvline(x=pas_alpha * (1 + redshift), ymin=0.65, ymax=0.85, lw='2', ls='--', color='r')
        ax.axvline(x=h2_1_0_s1 * (1 + redshift), ymin=0.65, ymax=0.85, lw='2', ls='--', color='r')
        ax.axvline(x=h2_1_0_s2 * (1 + redshift), ymin=0.65, ymax=0.85, lw='2', ls='--', color='r')
        ax.axvline(x=h2_1_0_s3 * (1 + redshift), ymin=0.65, ymax=0.85, lw='2', ls='--', color='r')
        ax.axvline(x=br_gamma * (1 + redshift), ymin=0.65, ymax=0.85, lw='2', ls='--', color='r')
        ax.axvline(x=br_delta * (1 + redshift), ymin=0.65, ymax=0.85, lw='2', ls='--', color='r')

        #axis_to_data = ax.transAxes + ax.transData.inverted()
        #axis_to_data.transform((xp, yp))  
        # this will allow you to transform any (xp, yp) axes point to
        # data coordinates. I don't need these lines for now but htey're useful.

        ax.text(x=pas_alpha * (1 + redshift), y=650.0, s=r'$\mathrm{Pas\ \alpha}$', fontsize=8)
        ax.text(x=h2_1_0_s2 * (1 + redshift), y=650.0, s=r'$\mathrm{H_2\ 1-0\, S(2)}$', fontsize=8)
        ax.text(x=h2_1_0_s3 * (1 + redshift), y=650.0, s=r'$\mathrm{H_2\ 1-0\, S(3)}$', fontsize=8)

        if h2_1_0_s1 * (1 + redshift) < 24100:
            ax.text(x=h2_1_0_s1 * (1 + redshift), y=650.0, s=r'$\mathrm{H_2\ 1-0\, S(1)}$', fontsize=8)

        if br_gamma * (1 + redshift) < 24100:
            ax.text(x=br_gamma * (1 + redshift), y=610.0, s=r'$\mathrm{Br\ \gamma}$', fontsize=8)

        if br_delta * (1 + redshift) < 24100:
            ax.text(x=br_delta * (1 + redshift), y=610.0, s=r'$\mathrm{Br\ \delta}$', fontsize=8)

    #ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))
    #ax.set_yscale('log')

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    fig.savefig(work_dir + obj_name + '_' + slitpos + '_' + band + '_tellinterp.png', dpi=300, bbox_inches='tight')

    return None

if __name__ == '__main__':
    
    obj_name = 'xl53'
    slitpos = 'BA'
    redshift = 0.031508

    ext_dir = '/Volumes/Bhavins_backup/ipac/Palomar_data/test/'
    date = ''

    obj_filename = ext_dir + date + obj_name + '_'\
     + slitpos + '_tellinterp_dispcor.fits'

    # give it the filename which has the dispersion corrected spectra
    obj_spec = fits.open(obj_filename)

    jsplt = obj_spec[0].header["WAT2_001"].split('=')[-1].split(' ')
    hsplt = obj_spec[0].header["WAT2_002"].split('=')[-1].split(' ')
    ksplt = obj_spec[0].header["WAT2_003"].split('=')[-1].split(' ')

    jstart = float(jsplt[4])
    delta_j = float(jsplt[5])

    hstart = float(hsplt[4])
    delta_h = float(hsplt[5])

    kstart = float(ksplt[4])
    delta_k = float(ksplt[5])

    num = 2048

    jwav = [jstart+i*delta_j for i in range(2048)]
    jspec = obj_spec[0].data[0,0]

    hwav = [hstart+i*delta_h for i in range(2048)]
    hspec = obj_spec[0].data[0,1]

    kwav = [kstart+i*delta_k for i in range(2048)]
    kspec = obj_spec[0].data[0,2]

    plotspec(ext_dir, date, jwav, jspec, obj_name, 'j', slitpos, redshift)
    plotspec(ext_dir, date, hwav, hspec, obj_name, 'h', slitpos, redshift)
    plotspec(ext_dir, date, kwav, kspec, obj_name, 'k', slitpos, redshift)

    sys.exit(0)