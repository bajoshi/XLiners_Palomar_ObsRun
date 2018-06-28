from __future__ import division

from astropy.io import fits
import numpy as np
from pyraf import iraf

import string
import os
import sys
import shutil

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
xliners_dir = home + '/Desktop/ipac/XLiners_Palomar_ObsRun/'

sys.path.append(xliners_dir)
import fix_lacosmic_header as fx
import telluric_interp_stellar_lines as ti
import div_obj_telluric as dv
import coadd_and_smooth as cs

def make_flat_fielded_images(work_dir, master_flat, sci_list):

    # change iraf working directory
    iraf.cd(work_dir)

    # loop over all images in sci list and flat field 
    # each one. Pyraf (while scripting) does not seem 
    # to like being given lists (e.g. @...)
    fh = open(work_dir + sci_list, 'r')

    for fl in fh.read().splitlines():
        sci_name = fl
        sci_f_name = sci_name.replace('.fits', '_f.fits')
        iraf.imarith(sci_name, '/', master_flat, sci_f_name)

    fh.close()

    return None

def get_list(work_dir, listname):

    fh = open(work_dir + listname, 'r')

    all_im = []
    for fl in fh.read().splitlines():
        all_im.append(fl)

    fh.close()

    return all_im

def do_a_minus_b(work_dir, sci_list, prefix):

    # change iraf working directory
    iraf.cd(work_dir)

    # first make sure that there are an integer 
    # number of ABBA sequences in the list

    # put all images in a python list first 
    all_im = get_list(work_dir, sci_list)
    all_im = [i.replace('.fits', '_f.fits') for i in all_im]

    # loop over all images and do A-B
    a_list = []
    b_list = []

    if len(all_im) % 4 == 0:

        for i in np.arange(0,len(all_im),4):
            
            im_a1 = all_im[i]
            im_b1 = all_im[i+1]
            im_b2 = all_im[i+2]
            im_a2 = all_im[i+3]

            a1_num = str(int(im_a1.split(prefix)[-1].split('_')[0]))
            a2_num = str(int(im_a2.split(prefix)[-1].split('_')[0]))
            b1_num = str(int(im_b1.split(prefix)[-1].split('_')[0]))
            b2_num = str(int(im_b2.split(prefix)[-1].split('_')[0]))

            iraf.imarith(im_a1, '-', im_b1, prefix + a1_num + '.fits')
            iraf.imarith(im_a2, '-', im_b2, prefix + a2_num + '.fits')

            iraf.imarith(im_b1, '-', im_a1, prefix + b1_num + '.fits')
            iraf.imarith(im_b2, '-', im_a2, prefix + b2_num + '.fits')

            a_list.append(prefix + a1_num + '.fits')
            a_list.append(prefix + a2_num + '.fits')

            b_list.append(prefix + b1_num + '.fits')
            b_list.append(prefix + b2_num + '.fits')

    else:
        print "Please check the length of the list and make"
        print "sure that you are giving an integer number of ABBA sequences"
        print "in the correct order. Exiting."
        sys.exit(0)

    return a_list, b_list

def combine_images(work_dir, a_list, b_list, obj_name):

    # change iraf working directory
    iraf.cd(work_dir)

    iraf.imcombine(a_list, obj_name + '_AB_median.fits', combine='median')
    iraf.imcombine(b_list, obj_name + '_BA_median.fits', combine='median')

    return

def save_AB_lists(work_dir, obj_name, a_list, b_list):

    fa = open(work_dir + obj_name + '_AB.lis', 'wa')
    fb = open(work_dir + obj_name + '_BA.lis', 'wa')

    for i in range(len(a_list)):
        fa.write(str(a_list[i]) + '\n')
        fb.write(str(b_list[i]) + '\n')

    fa.close()
    fb.close()

    return None

def copy_from_raw(work_dir, raw_dir, sci_list):

    # copy all the science images over
    fh = open(work_dir + sci_list, 'r')

    for fl in fh.read().splitlines():
        if not os.path.isfile(work_dir + fl):
            shutil.copy(raw_dir + fl, work_dir + fl)

    fh.close()

    # also copy la_cosmic
    if not os.path.isfile(work_dir + 'la_cosmic.pro'):
        shutil.copy(xliners_dir + 'la_cosmic.pro', work_dir + 'la_cosmic.pro')

    # also copy normalized master flat
    if not os.path.isfile(work_dir + 'master_flatN.fits'):
        shutil.copy(raw_dir + 'master_flatN.fits', work_dir + 'master_flatN.fits')

    return None

def finish_combine(work_dir, raw_dir, obj_name, prefix):

    # change iraf working directory
    iraf.cd(work_dir)

    # put in code to make master dark and flat
    
    # first copy sci images from raw directory
    sci_list = obj_name + '_sci.lis'
    copy_from_raw(work_dir, raw_dir, sci_list)

    # make flat fielded images
    flat_norm = 'master_flatN.fits'
    make_flat_fielded_images(work_dir, flat_norm, sci_list)

    # A-B subtraction
    a_list, b_list = do_a_minus_b(work_dir, sci_list, prefix)

    # save the A and B lists to a file in case you need them later 
    save_AB_lists(work_dir, obj_name, a_list, b_list)

    # Combine all A's and B's
    # use these following two lines if you have the A 
    # and B lists and need to proceed from combining images
    # a_list = get_list(work_dir, obj_name + '_AB.lis')
    # b_list = get_list(work_dir, obj_name + '_BA.lis')

    a_list = string.join(a_list, ',')
    b_list = string.join(b_list, ',')
    combine_images(work_dir, a_list, b_list, obj_name)

    return None

def run_fix_lacosmic_header_setjd(work_dir, obj_name, refspecA, refspecB, setjd_ref=False):
 
    # change iraf working directory
    iraf.cd(work_dir)

    # fix the header for hte median combined files
    filenameA = obj_name + '_AB_median'
    filenameB = obj_name + '_BA_median'

    fx.fix_la_cosmic_header(work_dir, filenameA)
    fx.fix_la_cosmic_header(work_dir, filenameB)

    # Now run setjd on the median combined files and the flat fields of the ref specs
    filenameA += '_cleaned.fits'
    filenameB += '_cleaned.fits'
    iraf.setjd(filenameA, date='DATE', epoch='')
    iraf.setjd(filenameB, date='DATE', epoch='')

    # this only needs to be done once
    # i.e. the first time in a night
    # probalby bet to do it in hte IRAF CL
    if setjd_ref:
        ffA = refspecA.replace('.fits', '_f.fits')
        ffB = refspecB.replace('.fits', '_f.fits')
        iraf.setjd(ffA, date='DATE', epoch='')
        iraf.setjd(ffB, date='DATE', epoch='')

        # now fix the headers for the ref specs
        fx.fix_la_cosmic_header(work_dir, refspecA.split('.')[0], for_setjd=True)
        fx.fix_la_cosmic_header(work_dir, refspecB.split('.')[0], for_setjd=True)

    return None

def finish_final(work_dir, obj_name, tell_name, refspecA, refspecB):

    # change iraf working directory
    iraf.cd(work_dir)

    # check if the telluric stellar abs features have already been 
    # interpolated over. If not then do the interpolation and save.
    tell_interp_A = tell_name + '_AB_' + 'stellarlines_interp_norm_bksub.fits'
    tell_interp_B = tell_name + '_BA_' + 'stellarlines_interp_norm_bksub.fits'

    if not os.path.isfile(work_dir + tell_interp_A):
        ti.do_all(work_dir, tell_name, 'AB')

    if not os.path.isfile(work_dir + tell_interp_B):
        ti.do_all(work_dir, tell_name, 'BA')

    # div object by the interpolated telluric spectrum to normalize
    # and make sure that the REFSPEC1 header keyword is set correctly
    dv.div_obj(work_dir, obj_name, tell_name, 'AB', refspecA)
    dv.div_obj(work_dir, obj_name, tell_name, 'BA', refspecB)

    return None

if __name__ == '__main__':

    # load all required packages
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.echelle(_doprint=0)
    iraf.onedspec(_doprint=0)

    # definitions
    work_dir = '/Volumes/Bhavins_backup/ipac/Palomar_data/2016/2016B/baj_work_night2/'
    raw_dir = '/Volumes/Bhavins_backup/ipac/Palomar_data/2016/2016B/2016OCT22/'
    obj_name = 'hd216219'
    redshift = 0.0
    telluric = 'hip5164'
    prefix = 'tspec'

    refspecA = 'tspec0054.fits'
    refspecB = 'tspec0055.fits'

    #finish_combine(work_dir, raw_dir, telluric, prefix)
    #finish_combine(work_dir, raw_dir, obj_name, prefix)
    #sys.exit(0)

    """
        Make sure you run LA_COSMIC in IDL before this next step is run.
        You need to run LA_COSMIC on the median combined images and on
        the reference images for wavelength calibration
        for the reference images this only has to be done for a single 
        A and a B file per night because the .ec files for these 
        reference spectra will simply be copied for other objects.
    """

    #run_fix_lacosmic_header_setjd(work_dir, obj_name, refspecA, refspecB)
    #sys.exit(0)
    
    """
        the first time you reduce data for a night you'll have to run doecslit
        see the reduction notes for that
        If you are not running this script for first time in a night then you don't 
        need doecslit. You can proceed directly to running apall inside iraf
        after running the above two functions for the sci list of a given object.
        make sure that the apall settings are as shown in the notes.
        Check the notes carefully for apall.
        MAKE SURE to run apall on the telluric as well.
    """

    finish_final(work_dir, obj_name, telluric, refspecA, refspecB)
    sys.exit(0)

    """
        RUN DISPCOR in IRAF now!
        Check reduction notes.
    """
    cs.stack_and_finish(work_dir, obj_name, redshift, smooth_width=1.0)
    sys.exit(0)