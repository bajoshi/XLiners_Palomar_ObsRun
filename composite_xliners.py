from __future__ import division

from astropy.io import fits
import numpy as np
from astropy.modeling import models, fitting

import os
import sys

import matplotlib.pyplot as plt

home = os.getenv('HOME')  # Does not have a trailing slash at the end
xliners_dir = home + '/Desktop/ipac/XLiners_Palomar_ObsRun/'

sys.path.append(xliners_dir)
import coadd_and_smooth as cs
import plot_final_spectra as pfs

def add_spec(lam_em, flux_em, flux_grid, lam_grid):

    for i in range(len(lam_grid)):

        if i < len(lam_grid)-1:
            diff = (lam_grid[i+1] - lam_grid[i]) / 2
        
        lam_idx = np.where((lam_em >= lam_grid[i]-diff) & (lam_em < lam_grid[i]+diff))[0]
        flux_grid[i].append(flux_em[lam_idx])

    return flux_grid

def makefig():

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('$\lambda\ [\mu m]$')
    ax.set_ylabel('$f_{\lambda}\ [\mathrm{W/m^{-2}/\mu m]}$')

    return fig, ax

def set_minorticks(ax):

    ax.minorticks_on()
    ax.tick_params('both', width=1, length=3, which='minor')
    ax.tick_params('both', width=1, length=4.7, which='major')

    return ax

def mark_emission_lines_kband(kflux_high, ax_k):

    # In microns
    pas_alpha = 1.8751  # air
    h2_1_0_s3 = 1.9576
    h2_1_0_s2 = 2.0338
    h2_1_0_s1 = 2.1218
    br_gamma = 2.1655  # air
    br_delta = 1.9446  # air

    ax_k.axvline(x=pas_alpha, ymin=0.35, ymax=0.65, lw='1.5', ls='--', color='r')
    ax_k.axvline(x=h2_1_0_s1, ymin=0.35, ymax=0.65, lw='1.5', ls='--', color='r')
    ax_k.axvline(x=h2_1_0_s2, ymin=0.35, ymax=0.65, lw='1.5', ls='--', color='r')
    ax_k.axvline(x=h2_1_0_s3, ymin=0.35, ymax=0.65, lw='1.5', ls='--', color='r')
    ax_k.axvline(x=br_gamma,  ymin=0.35, ymax=0.65, lw='1.5', ls='--', color='r')
    ax_k.axvline(x=br_delta,  ymin=0.35, ymax=0.65, lw='1.5', ls='--', color='r')

    ax_k_ypos1 = kflux_high * 0.85
    ax_k_ypos2 = kflux_high * 0.75

    ax_k.text(x=pas_alpha, y=ax_k_ypos1, s=r'$\mathrm{Pas\ \alpha}$', fontsize=11)
    ax_k.text(x=h2_1_0_s2, y=ax_k_ypos1, s=r'$\mathrm{H_2\ 1-0\, S(2)}$', fontsize=11)
    ax_k.text(x=h2_1_0_s3, y=ax_k_ypos2, s=r'$\mathrm{H_2\ 1-0\, S(3)}$', fontsize=11)

    if h2_1_0_s1 < k_high:
        ax_k.text(x=h2_1_0_s1, y=ax_k_ypos2, s=r'$\mathrm{H_2\ 1-0\, S(1)}$', fontsize=11)

    if br_gamma < k_high:
        ax_k.text(x=br_gamma, y=ax_k_ypos1, s=r'$\mathrm{Br\ \gamma}$', fontsize=11)

    if br_delta < k_high:
        ax_k.text(x=br_delta, y=ax_k_ypos1, s=r'$\mathrm{Br\ \delta}$', fontsize=11)

    return ax_k

def mark_emission_lines_jband(jflux_high, ax_j):

    # wav in microns
    pas_beta = 1.2822
    feii = 1.2570

    ax_j.axvline(x=pas_beta, ymin=0.35, ymax=0.65, lw='2', ls='--', color='r')
    ax_j.axvline(x=feii, ymin=0.35, ymax=0.65, lw='2', ls='--', color='r')

    ax_j_ypos1 = jflux_high * 0.65
    ax_j_ypos2 = jflux_high * 0.45

    ax_j.text(x=pas_beta, y=ax_j_ypos1, s=r'$\mathrm{Pas\ \beta}$', fontsize=11)
    ax_j.text(x=feii, y=ax_j_ypos1, s=r'$\mathrm{[FeII]}$', fontsize=11)

    return ax_j

def mark_emission_lines_hband(hflux_high, ax_h):

    feii = 1.6436  # wav in microns

    ax_h.axvline(x=feii, ymin=0.35, ymax=0.65, lw='2', ls='--', color='r')

    ax_h_ypos1 = hflux_high * 0.65
    ax_h_ypos2 = hflux_high * 0.45

    ax_h.text(x=feii, y=ax_h_ypos1, s=r'$\mathrm{[FeII]}$', fontsize=11)

    return ax_h

def plotspec(dirname_list, fname_list, redshift_list, lam_grid_j, lam_grid_h, lam_grid_k, \
    flux_grid_j, flux_grid_h, flux_grid_k, method):

    # make figures
    fig_j, ax_j = makefig()
    fig_h, ax_h = makefig()
    fig_k, ax_k = makefig()

    # Draw horizontal zero line
    ax_j.axhline(y=0,linestyle='--')
    ax_h.axhline(y=0,linestyle='--')
    ax_k.axhline(y=0,linestyle='--')

    # Smoothing
    #flux_grid_j = cs.smoothspectra(flux_grid_j, width=3.0)
    #flux_grid_h = cs.smoothspectra(flux_grid_h, width=3.0)
    #flux_grid_k = cs.smoothspectra(flux_grid_k, width=3.0)

    # plotting
    ax_j.plot(lam_grid_j, flux_grid_j, '-', color='k', linewidth=1, markeredgecolor='k', markersize=4, zorder=10)
    ax_h.plot(lam_grid_h, flux_grid_h, '-', color='k', linewidth=1, markeredgecolor='k', markersize=4, zorder=10)
    ax_k.plot(lam_grid_k, flux_grid_k, '-', color='k', linewidth=1, markeredgecolor='k', markersize=4, zorder=10)

    """
    # For plotting individual lines
    count = 0
    for i in range(len(dirname_list)):

        dirname = dirname_list[i]
        fname = fname_list[i]
        redshift = redshift_list[i]

        if count == 0:
            if continuum_flag == 'rescale':
                j_medarr, j_medval, h_medarr, h_medval, k_medarr, k_medval = rescale(dirname_list, fname_list, redshift_list)

        # prepare spectra based on user given choice
        if continuum_flag == 'contsub':
            j_em, h_em, k_em, jflux_em, hflux_em, kflux_em = prep_spectra(dirname, fname, redshift, subtract_continuum=True)
        elif continuum_flag == 'rescale':
            j_em, h_em, k_em, jflux_em, hflux_em, kflux_em = prep_spectra(dirname, fname, redshift, subtract_continuum=False)

            # rescale the spectra
            jflux_em = (jflux_em / j_medarr[u]) * j_medval
            hflux_em = (hflux_em / h_medarr[u]) * h_medval
            kflux_em = (kflux_em / k_medarr[u]) * k_medval

        plot_indiv = False
        if plot_indiv:
            # smooth indiv spectra only for plotting purposes
            jflux_em = cs.smoothspectra(jflux_em, width=5.0)
            hflux_em = cs.smoothspectra(hflux_em, width=5.0)
            kflux_em = cs.smoothspectra(kflux_em, width=5.0)

            ax_j.plot(j_em, jflux_em, ls='-', linewidth=1, color='lightgray')
            ax_h.plot(h_em, hflux_em, ls='-', linewidth=1, color='lightgray')
            ax_k.plot(k_em, kflux_em, ls='-', linewidth=1, color='lightgray')

        count += 1
    """

    ax_j = set_minorticks(ax_j)
    ax_h = set_minorticks(ax_h)
    ax_k = set_minorticks(ax_k)

    j_low, j_high, h_low, h_high, k_low, k_high = get_lims(restframe=True)

    ax_j.set_xlim(j_low, j_high)
    ax_h.set_xlim(h_low, h_high)
    ax_k.set_xlim(k_low, k_high)

    # Get limits for plotting
    if 'median' in method:
        jflux_low, jflux_high = pfs.get_flux_lims(lam_grid_j, flux_grid_j, j_low, j_high, force_zero=False, extension_factor=0.1)
        hflux_low, hflux_high = pfs.get_flux_lims(lam_grid_h, flux_grid_h, h_low, h_high, force_zero=False, extension_factor=0.1)
        kflux_low, kflux_high = pfs.get_flux_lims(lam_grid_k, flux_grid_k, k_low, k_high, force_zero=False, extension_factor=0.1)
    elif 'mean' in method:
        jflux_low, jflux_high = -30, 40
        hflux_low, hflux_high = -60, 70
        kflux_low, kflux_high = -50, 80

    ax_j.set_ylim(jflux_low, jflux_high)
    ax_h.set_ylim(hflux_low, hflux_high)
    ax_k.set_ylim(kflux_low, kflux_high)

    # add vertical lines to show emission lines
    ax_j = mark_emission_lines_jband(jflux_high, ax_j)
    ax_h = mark_emission_lines_hband(hflux_high, ax_h)
    ax_k = mark_emission_lines_kband(kflux_high, ax_k)

    # save figure
    fig_j.savefig(home + '/Desktop/ipac/composite_and_indiv_spec_plots/' \
        + 'j_fullstack_' + method + '.eps', dpi=300, bbox_inches='tight')
    fig_h.savefig(home + '/Desktop/ipac/composite_and_indiv_spec_plots/' \
        + 'h_fullstack_' + method + '.eps', dpi=300, bbox_inches='tight')
    fig_k.savefig(home + '/Desktop/ipac/composite_and_indiv_spec_plots/' \
        + 'k_fullstack_' + method + '.eps', dpi=300, bbox_inches='tight')

    # Save the stacked spectra
    np.savetxt(home + '/Desktop/ipac/composite_xliners_spectra/' + 'j_fullstack.txt', \
        np.c_[lam_grid_j, flux_grid_j], fmt='%.4e', delimiter=' ', header='wav[micron] spec[W/m2/micron]')
    np.savetxt(home + '/Desktop/ipac/composite_xliners_spectra/' + 'h_fullstack.txt', \
        np.c_[lam_grid_h, flux_grid_h], fmt='%.4e', delimiter=' ', header='wav[micron] spec[W/m2/micron]')
    np.savetxt(home + '/Desktop/ipac/composite_xliners_spectra/' + 'k_fullstack.txt', \
        np.c_[lam_grid_k, flux_grid_k], fmt='%.4e', delimiter=' ', header='wav[micron] spec[W/m2/micron]')

    #plt.show()

    return None

def rescale(dirname_list, fname_list, redshift_list):

    j_medarr = np.zeros(len(dirname_list))
    h_medarr = np.zeros(len(dirname_list))
    k_medarr = np.zeros(len(dirname_list))

    for i in range(len(dirname_list)):

        dirname = dirname_list[i]
        fname = fname_list[i]
        redshift = redshift_list[i]

        j_em, h_em, k_em, jflux_em, hflux_em, kflux_em = prep_spectra(dirname, fname, redshift, subtract_continuum=False)

        j_idx = np.where((j_em >= 10000) & (j_em <= 11000))[0]
        h_idx = np.where((h_em >= 14000) & (h_em <= 15000))[0]
        k_idx = np.where((k_em >= 19000) & (j_em <= 20000))[0]

        j_medarr[i] = np.median(jflux_em[j_idx])
        h_medarr[i] = np.median(hflux_em[h_idx])
        k_medarr[i] = np.median(kflux_em[k_idx])

    j_medval = np.median(j_medarr)
    h_medval = np.median(h_medarr)
    k_medval = np.median(k_medarr)

    return j_medarr, j_medval, h_medarr, h_medval, k_medarr, k_medval

def prep_spectra(dirname, obj_name, redshift, subtract_continuum=True):

    # read in spectra
    filename_j = dirname + '/' + obj_name + '_Jspec_phys_units.txt'
    filename_h = dirname + '/' + obj_name + '_Hspec_phys_units.txt'
    filename_k = dirname + '/' + obj_name + '_Kspec_phys_units.txt'

    specfile_j = np.genfromtxt(filename_j, dtype=None, names=True)
    specfile_h = np.genfromtxt(filename_h, dtype=None, names=True)
    specfile_k = np.genfromtxt(filename_k, dtype=None, names=True)

    jwav = specfile_j['wavmicron']
    jspec = specfile_j['specWm2micron']

    hwav = specfile_h['wavmicron']
    hspec = specfile_h['specWm2micron']

    kwav = specfile_k['wavmicron']
    kspec = specfile_k['specWm2micron']

    # Subtract continuum
    if subtract_continuum:
        # fit polynomial and subtract continuum
        # make sure that the average across stacking 
        # wavelength range is close to zero after continuum subtraction

        # initialize polynomail and fit
        p_init_jh = models.Polynomial1D(degree=3)
        p_init_k = models.Polynomial1D(degree=3)
        fit_p = fitting.LinearLSQFitter()

        # only fit polynomial to range with believalbe fluxes
        j_low, j_high, h_low, h_high, k_low, k_high = get_lims()
        # It seems like some of my reductions aren't correct and still have 
        # atmospheric lines in them. I'll have to redo the redcutions for these objects.

        # Find indices and fit
        j_low_idx = np.argmin(abs(jwav - j_low))
        j_high_idx = np.argmin(abs(jwav - j_high))

        h_low_idx = np.argmin(abs(hwav - h_low))
        h_high_idx = np.argmin(abs(hwav - h_high))

        k_low_idx = np.argmin(abs(kwav - k_low))
        k_high_idx = np.argmin(abs(kwav - k_high))

        pj = fit_p(p_init_jh, jwav[j_low_idx:j_high_idx], jspec[j_low_idx:j_high_idx])
        ph = fit_p(p_init_jh, hwav[h_low_idx:h_high_idx], hspec[h_low_idx:h_high_idx])
        pk = fit_p(p_init_k, kwav[k_low_idx:k_high_idx], kspec[k_low_idx:k_high_idx])

        # Check fit
        #check_polynomial_fit(pj, jwav, jspec)
        #check_polynomial_fit(ph, hwav, hspec)
        #check_polynomial_fit(pk, kwav, kspec)

        # subtract continuum
        jspec -= pj(jwav)
        hspec -= ph(hwav)
        kspec -= pk(kwav)

        # check mean after subtracting continuum
        # this should be close to zero
        #print np.mean(jspec[j_low_idx:j_high_idx])
        #print np.mean(hspec[h_low_idx:h_high_idx])
        #print np.mean(kspec[k_low_idx:k_high_idx])

    # deredshift each spectrum and put in array
    j_em = jwav / (1 + redshift)
    h_em = hwav / (1 + redshift)
    k_em = kwav / (1 + redshift)

    # If your fluxes are correct then the flux will also have to be deredshifted
    # right now these lines for flux deredshifting don't matter
    jflux_em = jspec * (1 + redshift)
    hflux_em = hspec * (1 + redshift)
    kflux_em = kspec * (1 + redshift)

    return j_em, h_em, k_em, jflux_em, hflux_em, kflux_em

def check_polynomial_fit(poly, wav, spec):

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(wav, spec, color='k', lw=2)
    ax.plot(wav, poly(wav), color='b', lw=1)

    plt.show()

    return None

def get_lims(restframe=False):

    j_low = 11500
    j_high = 13500

    h_low = 14900
    h_high = 18000

    k_low = 19500
    k_high = 24100

    # uncomment if you want to see the entire grid
    # right now its only showing the region with important lines
    #if restframe:
    #    j_low = 9800
    #    j_high = 13200

    #    h_low = 12800
    #    h_high = 17400

    #    k_low = 16800
    #    k_high = 23400

    if restframe:
        j_low = 11000
        j_high = 13000

        h_low = 14200
        h_high = 17000

        k_low = 18600
        k_high = 23000

    # Convert to microns 
    j_low /= 1e4 
    j_high /= 1e4
    h_low /= 1e4
    h_high /= 1e4
    k_low /= 1e4
    k_high /= 1e4

    return j_low, j_high, h_low, h_high, k_low, k_high

if __name__ == '__main__':

    # read in user command line options
    # sys.argv[0] is the name of hte script
    try:
        continuum_flag = 'contsub'
        combine_flag = sys.argv[2]

    except IndexError as e:
        print e
        print "This code needs a method (rescaling all continua or"
        print "subtracting individual continua before combining) "
        print "to be set for the individual continua."
        print "Keyword options are -- rescale or contsub"
        print "It also needs a combining method -- median or mean."
        print "E.g. >> python composite_xliners.py rescale mean"
        sys.exit(0)

    # definitions
    obj_name_arr = np.array(['xl49','xl53','xl55','xl100','xl124',\
    'xl208_total','xl229','xl435','xl692','xw3','xw244','xw661', 'xw588', 'xw867'])
    z_arr = np.array([0.1647, 0.031508, 0.131864, 0.0821633, 0.0875584,\
    0.172875, 0.149458, 0.0848958, 0.0806426, 0.14795, 0.13154, 0.04733, 0.03701, 0.07083])

    ## for LOW Z stack
    low_z_stack = False
    low_z_lim = 0.2
    if low_z_stack:

        z_list = []
        obj_list = []
        for v in range(len(obj_name_arr)):
            
            if z_arr[v] <= low_z_lim:
                print z_arr[v], obj_name_arr[v]
                z_list.append(z_arr[v])
                obj_list.append(obj_name_arr[v])

        z_arr = np.array(z_list)
        obj_name_arr = np.array(obj_list)

    comp_dir = home + '/Desktop/ipac/composite_xliners_spectra/'

    # create lambda and flux grid for all bands
    # observed frame limits
    lowest_z = np.min(z_arr)
    highest_z = np.max(z_arr)
    j_low, j_high, h_low, h_high, k_low, k_high = get_lims()

    # rest frame limits
    # The endpoints are chosen based on hte limits of the observed 
    # band and the lowest and highest redshift of the sample.
    j_grid_low = 0.98  #np.floor(j_low / (1+highest_z))
    j_grid_high = 1.31  #np.ceil(j_high / (1+lowest_z))

    h_grid_low = 1.27  #np.floor(h_low / (1+highest_z))
    h_grid_high = 1.75  #np.ceil(h_high / (1+lowest_z))

    k_grid_low = 1.66  #np.floor(k_low / (1+highest_z))
    k_grid_high = 2.34  #np.ceil(k_high / (1+lowest_z))
    # I've commented out the floor and ceil lines because
    # they worked fine when the numbers were in angstroms
    # but now that the wavelenghts are in microns they 
    # don't work right. This is because these functions 
    # will always return an integer answer. 

    # rest frame grids
    # I'm using linspace because that will include the endpoints
    # I chose to have 1500 points by trial and error
    total_lam_grid_elem = 2086
    lam_grid_j = np.linspace(j_grid_low, j_grid_high, num=total_lam_grid_elem)
    lam_grid_h = np.linspace(h_grid_low, h_grid_high, num=total_lam_grid_elem)
    lam_grid_k = np.linspace(k_grid_low, k_grid_high, num=total_lam_grid_elem)

    flux_grid_j = np.empty(len(lam_grid_j))
    flux_grid_h = np.empty(len(lam_grid_h))
    flux_grid_k = np.empty(len(lam_grid_k))

    # convert to list to be able to set up each 
    # individual element as an empty list itself.
    # I need the numpy aray first to get an array 
    # of known length.
    flux_grid_j = flux_grid_j.tolist()
    flux_grid_h = flux_grid_h.tolist()
    flux_grid_k = flux_grid_k.tolist()

    # loop over all objects
    # first make a proper list of all object paths
    """
    dirname_list = []
    fname_list = []
    redshift_list = []
    for dirname, subdirname, flist in os.walk(comp_dir):
        for fname in flist:

            curr_obj = fname.split('_')[0]

            if curr_obj == 'xl208':
                curr_obj = 'xl208_total'

            if curr_obj in obj_name_arr:
                curr_z = z_arr[np.where(obj_name_arr == curr_obj)[0]]

                # skipping first dir because its the mac default .DS_Store directory
                if fname == '.DS_Store':
                    continue

                dirname_list.append(dirname)
                fname_list.append(fname)
                redshift_list.append(curr_z[0])
    """

    # Full list
    dates_list = ['20160523', '20160523', \
    '20160524', '20160524', '20160524', '20160524', '20160524', '20160524', '20160524', '20160524', \
    '20161021', '20161021', '20161021', '20161021', '20161021', '20161021', '20161021', '20161021', \
    '20161022', '20161022', '20161022', '20161022', '20161022', '20161022', '20161022', '20161022', '20161022', '20161022', \
    '20170509', '20170509', '20170509', '20170509', '20170509', '20170509',\
    '20170511', '20170511', '20170511', '20170511', '20170511', '20170511', '20170511', '20170511']
    obj_ondate_list = ['xl53_A', 'xl53_B', \
    'xl53_A', 'xl53_B', 'xl124_A', 'xl124_B', 'xl208_total_A', 'xl208_total_B', 'xl229_A', 'xl229_B', \
    'xw3_A', 'xw3_B', 'xw244_A', 'xw244_B', 'xw588_A', 'xw588_B', 'xw867_A', 'xw867_B', \
    'xw3_A', 'xw3_B', 'xw244_A', 'xw244_B', 'xw588_A', 'xw588_B', 'xw661_A', 'xw661_B', 'xw867_A', 'xw867_B', \
    'xl53_A', 'xl53_B', 'xl55_A', 'xl55_B', 'xl229_A', 'xl229_B', \
    'xl53_A', 'xl53_B', 'xl100_A', 'xl100_B', 'xl435_A', 'xl435_B', 'xl692_A', 'xl692_B']

    """ # small list
    dates_list = ['20160523', '20160523', \
    '20160524', '20160524', '20160524', '20160524', '20160524', '20160524', '20160524', '20160524', \
    '20161021', '20161021', '20161021', '20161021', \
    '20161022', '20161022', '20161022', '20161022', '20161022', '20161022', \
    '20170509', '20170509', '20170509', '20170509', '20170509', '20170509',\
    '20170511', '20170511', '20170511', '20170511', '20170511', '20170511', '20170511', '20170511']
    obj_ondate_list = ['xl53_A', 'xl53_B', \
    'xl53_A', 'xl53_B', 'xl124_A', 'xl124_B', 'xl208_total_A', 'xl208_total_B', 'xl229_A', 'xl229_B', \
    'xw3_A', 'xw3_B', 'xw244_A', 'xw244_B', \
    'xw3_A', 'xw3_B', 'xw244_A', 'xw244_B', 'xw661_A', 'xw661_B', \
    'xl53_A', 'xl53_B', 'xl55_A', 'xl55_B', 'xl229_A', 'xl229_B', \
    'xl53_A', 'xl53_B', 'xl100_A', 'xl100_B', 'xl435_A', 'xl435_B', 'xl692_A', 'xl692_B']
    """

    # create redshift list
    redshift_list = []
    dirname_list = []
    for w in range(len(obj_ondate_list)):
        curr_obj = obj_ondate_list[w].split('_')[0]
        if curr_obj == 'xl208':
            curr_obj = 'xl208_total'

        curr_z = z_arr[np.where(obj_name_arr == curr_obj)[0]]
        #print curr_obj, curr_z
        redshift_list.append(curr_z[0])
        dirname_list.append(home + '/Desktop/ipac/composite_xliners_spectra/' + dates_list[w])

    # now actually looping over objects to stack
    count = 0
    for u in range(len(obj_ondate_list)):

        dirname = dirname_list[u]
        fname = obj_ondate_list[u]
        curr_z = redshift_list[u]

        print "Currently at object:", fname, "at redshift:", curr_z

        # set up spectrum grid 
        # only needs to be done once 
        # --> on the first time the loop iterates.
        # you also only need to run the rescaling routine only once
        if count == 0:
            for i in range(total_lam_grid_elem):
                flux_grid_j[i] = []
                flux_grid_h[i] = []
                flux_grid_k[i] = []

            if continuum_flag == 'rescale':
                j_medarr, j_medval, h_medarr, h_medval, k_medarr, k_medval = rescale(dirname_list, fname_list, redshift_list)
        
        # prepare spectra based on user given choice
        if continuum_flag == 'contsub':
            j_em, h_em, k_em, jflux_em, hflux_em, kflux_em = prep_spectra(dirname, fname, curr_z, subtract_continuum=True)
        elif continuum_flag == 'rescale':
            j_em, h_em, k_em, jflux_em, hflux_em, kflux_em = prep_spectra(dirname, fname, curr_z, subtract_continuum=False)

            # rescale the spectra
            jflux_em = (jflux_em / j_medarr[u]) * j_medval
            hflux_em = (hflux_em / h_medarr[u]) * h_medval
            kflux_em = (kflux_em / k_medarr[u]) * k_medval

        # add spec to empty spec which is aligned with the lam grid
        # based on add_spec in grid_coadd.py
        flux_grid_j = add_spec(j_em, jflux_em, flux_grid_j, lam_grid_j)
        flux_grid_h = add_spec(h_em, hflux_em, flux_grid_h, lam_grid_h)
        flux_grid_k = add_spec(k_em, kflux_em, flux_grid_k, lam_grid_k)

        count += 1

    # take median or mean of all points appended from indiv spectra
    if combine_flag == 'median':

        for l in range(len(lam_grid_j)):

            flux_grid_j[l] = np.median(np.concatenate(flux_grid_j[l]).ravel())
            flux_grid_h[l] = np.median(np.concatenate(flux_grid_h[l]).ravel())
            flux_grid_k[l] = np.median(np.concatenate(flux_grid_k[l]).ravel())
    
    elif combine_flag == 'mean':

        for l in range(len(lam_grid_j)):

            flux_grid_j[l] = np.mean(np.concatenate(flux_grid_j[l]).ravel())
            flux_grid_h[l] = np.mean(np.concatenate(flux_grid_h[l]).ravel())
            flux_grid_k[l] = np.mean(np.concatenate(flux_grid_k[l]).ravel())

    # plot all spectra and stacks for all three bands
    if low_z_stack:
        plotspec(dirname_list, obj_ondate_list, redshift_list, lam_grid_j, lam_grid_h, lam_grid_k, \
            flux_grid_j, flux_grid_h, flux_grid_k, method=continuum_flag + '_' + combine_flag + '_' + str(low_z_lim).replace('.','p'))
    else:
        plotspec(dirname_list, obj_ondate_list, redshift_list, lam_grid_j, lam_grid_h, lam_grid_k, \
            flux_grid_j, flux_grid_h, flux_grid_k, method=continuum_flag + '_' + combine_flag)

    sys.exit(0)
