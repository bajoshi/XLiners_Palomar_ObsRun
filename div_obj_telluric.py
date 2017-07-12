from __future__ import division

from astropy.io import fits
import numpy as np

import os
import sys

import matplotlib.pyplot as plt 

if __name__ == '__main__':

    tell_name = 'hip75230'
    obj_name = 'xl55'
    slitpos = 'BA'

    ext_dir = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/'
    date = '20170509/'

    # read in object spec
    obj_filename = ext_dir + date + obj_name + '_' + slitpos + '_bksub.fits'
    obj_spec = fits.open(obj_filename)

    # read in telluric spectrum iwth the stellar lines interpolated over
    tell_filename = ext_dir + date + tell_name + '_' + slitpos + '_stellarlines_interp_norm_bksub.fits'
    tell_spec_hdu = fits.open(tell_filename)

    tell_jspec_norm = tell_spec_hdu[0].data[0,0]
    tell_hspec_norm = tell_spec_hdu[0].data[0,1]
    tell_kspec_norm = tell_spec_hdu[0].data[0,2]

    # read in object spectrum and divide by telluric standard star and write
    obj_spec[0].data[0,0] /= tell_jspec_norm
    obj_spec[0].data[0,1] /= tell_hspec_norm
    obj_spec[0].data[0,2] /= tell_kspec_norm

    try:
        obj_spec[0].header['REFSPEC1'] += '.ec'
    except KeyError as e:
        refspec = raw_input('What is the ref spec for the ' + slitpos + ' position for this object? ')
        obj_spec[0].header.set('REFSPEC1', refspec)

    obj_spec.writeto(obj_filename.replace('.fits', '_tellinterp_norm.fits'), clobber=True)

    sys.exit(0)