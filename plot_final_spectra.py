from __future__ import division

from astropy.io import fits
import numpy as np

import os
import sys

import matplotlib.pyplot as plt 

def plotspec(wav, spec, obj_name, band, slitpos):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$\mathrm{\lambda\ [\AA]}$')
    ax.set_ylabel(r'$\mathrm{f\ [arbitrary scale]}$')

    ax.plot(wav, spec, color='k')

    if band == 'j':
        ax.set_xlim(11500, 13500)
        ax.set_ylim(0, 50)
    elif band == 'h':
        ax.set_xlim(14900, 18000)
        ax.set_ylim(0, 350)
    elif band == 'k':
        ax.set_xlim(19600, 24100)
        ax.set_ylim(0, 350)

    #ax.get_yaxis().get_major_formatter().set_powerlimits((0, 0))
    #ax.set_yscale('log')

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')
    ax.grid(True)

    obsdir = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/20170511/'

    fig.savefig(obsdir + obj_name + '_' + slitpos + '_' + band + '_tellinterp.png', dpi=300, bbox_inches='tight')

    return None

if __name__ == '__main__':
    
    obj_name = 'xl692'
    slitpos = 'BA'
    obj_filename = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/20170511/' + obj_name + '_'\
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

    plotspec(jwav, jspec, obj_name, 'j', slitpos)
    plotspec(hwav, hspec, obj_name, 'h', slitpos)
    plotspec(kwav, kspec, obj_name, 'k', slitpos)

    sys.exit(0)