from __future__ import division

from astropy.io import fits
import numpy as np
from scipy.interpolate import griddata

import os
import sys

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
xliners_dir = home + '/Desktop/ipac/XLiners_Palomar_ObsRun/'

if __name__ == '__main__':

    # set directories and names
    palomar_datadir = '/Volumes/Bhavins_backup/ipac/Palomar_data/'
    workdir = '2016/2016A/P2016A/fluxcal_workdir2018_2016Anight2/'
    obj_name = 'hd216219'

    # Read in A and B spectra for J,H,K stored in plain text files
    ja = np.genfromtxt(palomar_datadir + workdir + obj_name + '_A_Jspec_dispcor.txt', dtype=None, names=True)
    ha = np.genfromtxt(palomar_datadir + workdir + obj_name + '_A_Hspec_dispcor.txt', dtype=None, names=True)
    ka = np.genfromtxt(palomar_datadir + workdir + obj_name + '_A_Kspec_dispcor.txt', dtype=None, names=True)

    jb = np.genfromtxt(palomar_datadir + workdir + obj_name + '_B_Jspec_dispcor.txt', dtype=None, names=True)
    hb = np.genfromtxt(palomar_datadir + workdir + obj_name + '_B_Hspec_dispcor.txt', dtype=None, names=True)
    kb = np.genfromtxt(palomar_datadir + workdir + obj_name + '_B_Kspec_dispcor.txt', dtype=None, names=True)

    # save arrays
    ja_wav = ja['wavA']
    ha_wav = ha['wavA']
    ka_wav = ka['wavA']
    ja_spec = ja['specDNs']
    ha_spec = ha['specDNs']
    ka_spec = ka['specDNs']

    jb_wav = jb['wavA']
    hb_wav = hb['wavA']
    kb_wav = kb['wavA']
    jb_spec = jb['specDNs']
    hb_spec = hb['specDNs']
    kb_spec = kb['specDNs']

    # Wavelength is in Angstrom. Convert to micron
    ja_wav /= 1e4
    ha_wav /= 1e4
    ka_wav /= 1e4

    jb_wav /= 1e4
    hb_wav /= 1e4
    kb_wav /= 1e4

    # create interpolation wavelength grids
    #jgrid = np.arange()
    #hgrid = np.arange()

    kstep = 0.0002875  # in microns
    samp = 2086
    kstart = 1.88
    kgrid = np.linspace(start=kstart, stop=kstart + samp*kstep, num=samp, endpoint=False)

    # Now interpolate
    ka_spec_interp = griddata(points=ka_wav, values=ka_spec, xi=kgrid, method='linear')

    # Plot to check
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(kgrid, ka_spec_interp)
    ax.plot(ka_wav, ka_spec+10, color='k')

    plt.show()
    """

    # Divide by RSRF (Relative Spectral Response Function)
    # Read in RSRF
    k_rsrf = np.genfromtxt()

    sys.exit(0)