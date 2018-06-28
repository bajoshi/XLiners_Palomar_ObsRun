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
    workdir = '2016/2016B/baj_work_night1/'
    obj_name = 'hd216219'

    # Read in RSRF
    j_rsrf = np.genfromtxt(palomar_datadir + '2016/2016A/P2016A/fluxcal_workdir2018_2016Anight2/' + 'RSRF_HD216219jnew.txt', \
        dtype=None, names=['wav', 'rsrf'])
    k_rsrf = np.genfromtxt(palomar_datadir + '2016/2016A/P2016A/fluxcal_workdir2018_2016Anight2/' + 'RSRF_HD216219knew.txt', \
        dtype=None, names=['wav', 'rsrf'])

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

    # Plot spectra to check
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(jb_wav, jb_spec)
    ax.axhline(y=0.0, color='k', ls='--')
    ax.minorticks_on()
    plt.show()
    sys.exit(0)
    """

    # Wavelength is in Angstrom. Convert to micron
    ja_wav /= 1e4
    ha_wav /= 1e4
    ka_wav /= 1e4

    jb_wav /= 1e4
    hb_wav /= 1e4
    kb_wav /= 1e4

    # create interpolation wavelength grids
    jstep = 0.000172784 # in microns
    samp = 2086
    jstart = 1.12955
    #jgrid = np.linspace(start=jstart, stop=jstart + samp*jstep, num=samp, endpoint=False)

    jgrid = j_rsrf['wav']

    kstep = 0.0002875  # in microns
    kstart = 1.88
    kgrid = np.linspace(start=kstart, stop=kstart + samp*kstep, num=samp, endpoint=False)

    # Now interpolate
    # Both A and B spectra are interpolated to the same grid
    ja_spec_interp = griddata(points=ja_wav, values=ja_spec, xi=jgrid, method='linear')
    jb_spec_interp = griddata(points=jb_wav, values=jb_spec, xi=jgrid, method='linear')

    ka_spec_interp = griddata(points=ka_wav, values=ka_spec, xi=kgrid, method='linear')
    kb_spec_interp = griddata(points=kb_wav, values=kb_spec, xi=kgrid, method='linear')

    # Plot to check
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(kgrid, ka_spec_interp)
    ax.plot(ka_wav, ka_spec+10, color='k')
    plt.show()
    """

    # Divide by RSRF (Relative Spectral Response Function)
    # Divde new grided spectrum by RSRF to get spectrum in physical units # W/m^2/miron
    ja_spec_interp /= j_rsrf['rsrf']
    jb_spec_interp /= j_rsrf['rsrf']

    ka_spec_interp /= k_rsrf['rsrf']
    kb_spec_interp /= k_rsrf['rsrf']

    # save new physical spectrum
    # Jband
    np.savetxt(palomar_datadir + workdir + obj_name + '_A_Jspec_phys_units.txt', \
        np.c_[jgrid, ja_spec_interp], fmt='%.4e', delimiter=' ', header='wav[micron] spec[W/m2/micron]')
    np.savetxt(palomar_datadir + workdir + obj_name + '_B_Jspec_phys_units.txt', \
        np.c_[jgrid, jb_spec_interp], fmt='%.4e', delimiter=' ', header='wav[micron] spec[W/m2/micron]')

    # Kband
    np.savetxt(palomar_datadir + workdir + obj_name + '_A_Kspec_phys_units.txt', \
        np.c_[kgrid, ka_spec_interp], fmt='%.4e', delimiter=' ', header='wav[micron] spec[W/m2/micron]')
    np.savetxt(palomar_datadir + workdir + obj_name + '_B_Kspec_phys_units.txt', \
        np.c_[kgrid, kb_spec_interp], fmt='%.4e', delimiter=' ', header='wav[micron] spec[W/m2/micron]')

    # plot new spectrum
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(kgrid, ka_spec_interp-3e-15)
    ax.plot(kgrid, kb_spec_interp)
    plt.show()

    sys.exit(0)