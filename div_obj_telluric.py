from __future__ import division

from astropy.io import fits
import numpy as np

import os
import sys

import matplotlib.pyplot as plt 

def div_obj(work_dir, obj_name, tell_name, slitpos, refspec):

    # read in object spec
    obj_filename = work_dir + obj_name + '_' + slitpos + '_bksub.fits'
    obj_spec = fits.open(obj_filename)

    # read in telluric spectrum iwth the stellar lines interpolated over
    tell_filename = work_dir + tell_name + '_' + slitpos + '_stellarlines_interp_norm_bksub.fits'
    tell_spec_hdu = fits.open(tell_filename)

    tell_jspec_norm = tell_spec_hdu[0].data[0,0]
    tell_hspec_norm = tell_spec_hdu[0].data[0,1]
    tell_kspec_norm = tell_spec_hdu[0].data[0,2]

    # read in object spectrum and divide by telluric standard star and write
    obj_spec[0].data[0,0] /= tell_jspec_norm
    obj_spec[0].data[0,1] /= tell_hspec_norm
    obj_spec[0].data[0,2] /= tell_kspec_norm

    # set refspec explicity
    if len(refspec.split('.')) > 1:

        if 'cleaned' not in refspec:
            refspec = refspec.replace('.fits', '_cleaned.fits')

        refspec = refspec.replace('.fits', '.ec')

    elif len(refspec.split['.']) == 1:
        refspec += '_cleaned.ec'

    obj_spec[0].header.set('REFSPEC1', refspec)

    obj_spec.writeto(obj_filename.replace('.fits', '_tellinterp_norm.fits'), overwrite=True)

    return None

if __name__ == '__main__':

    tell_name = 'hd203856'
    obj_name = 'sqas1'
    slitpos = 'AB'
    refspec = 'SQ0043.fits'

    ext_dir = '/Volumes/Bhavins_backup/ipac/Palomar_data/2009/work/'

    div_obj(ext_dir, obj_name, tell_name, slitpos, refspec)

    sys.exit(0)