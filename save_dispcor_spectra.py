from __future__ import division

from astropy.io import fits
import numpy as np

import os
import sys

home = os.getenv('HOME')  # Does not have a trailing slash at the end
xliners_dir = home + '/Desktop/ipac/XLiners_Palomar_ObsRun/'

sys.path.append(xliners_dir)
import coadd_and_smooth as cs

if __name__ == '__main__':
    
    # set directories and names
    palomar_datadir = '/Volumes/Bhavins_backup/ipac/Palomar_data/'
    workdir = '2016/2016B/baj_work_night2/'
    obj_name = 'hd216219'
    exptime = 120.0

    # open discor given spectrum 
    hdu_A = fits.open(palomar_datadir + workdir + obj_name + '_AB_tellinterp_dispcor.fits')
    hdu_B = fits.open(palomar_datadir + workdir + obj_name + '_BA_tellinterp_dispcor.fits')

    # Now get the spectra
    jspec_ab, hspec_ab, kspec_ab, jwav_ab, hwav_ab, kwav_ab, jstart_ab, hstart_ab, kstart_ab = cs.get_spectra(hdu_A)
    jspec_ba, hspec_ba, kspec_ba, jwav_ba, hwav_ba, kwav_ba, jstart_ba, hstart_ba, kstart_ba = cs.get_spectra(hdu_B)

    # Divide by exposure time to get units to DN/s
    jspec_ab /= exptime
    hspec_ab /= exptime
    kspec_ab /= exptime

    jspec_ba /= exptime
    hspec_ba /= exptime
    kspec_ba /= exptime

    # Save the spectra as plain text files
    # All A spectra
    np.savetxt(palomar_datadir + workdir + obj_name + '_A_Jspec_dispcor.txt', \
        np.c_[jwav_ab, jspec_ab], fmt='%.4e', delimiter=' ', header='wav[A] spec[DN/s]')
    np.savetxt(palomar_datadir + workdir + obj_name + '_A_Hspec_dispcor.txt', \
        np.c_[hwav_ab, hspec_ab], fmt='%.4e', delimiter=' ', header='wav[A] spec[DN/s]')
    np.savetxt(palomar_datadir + workdir + obj_name + '_A_Kspec_dispcor.txt', \
        np.c_[kwav_ab, kspec_ab], fmt='%.4e', delimiter=' ', header='wav[A] spec[DN/s]')

    # All B spectra
    np.savetxt(palomar_datadir + workdir + obj_name + '_B_Jspec_dispcor.txt', \
        np.c_[jwav_ba, jspec_ba], fmt='%.4e', delimiter=' ', header='wav[A] spec[DN/s]')
    np.savetxt(palomar_datadir + workdir + obj_name + '_B_Hspec_dispcor.txt', \
        np.c_[hwav_ba, hspec_ba], fmt='%.4e', delimiter=' ', header='wav[A] spec[DN/s]')
    np.savetxt(palomar_datadir + workdir + obj_name + '_B_Kspec_dispcor.txt', \
        np.c_[kwav_ba, kspec_ba], fmt='%.4e', delimiter=' ', header='wav[A] spec[DN/s]')

    # close HDUs
    hdu_A.close()
    hdu_B.close()

    sys.exit(0)