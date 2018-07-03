from __future__ import division

import numpy as np

import os
import sys
import time
import datetime

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, AnchoredText

home = os.getenv('HOME')  # Does not have a trailing slash at the end
xliners_dir = home + '/Desktop/ipac/XLiners_Palomar_ObsRun/'

sys.path.append(xliners_dir)
import bpt_plots as bpt

if __name__ == '__main__':
    
    # Read in data
    cat = np.genfromtxt(taffydir + 'targets_small', dtype=None, delimiter=',', names=True)

    oiii_hbeta = cat['lgoiii2hb']
    oi_halpha = cat['lgoi2ha']

    # Plot 
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_ylabel(r'$\mathrm{log\left( \frac{[OIII]}{H\beta} \right)}$', fontsize=15)
    ax.set_xlabel(r'$\mathrm{log\left( \frac{[OI]}{H\alpha} \right)}$', fontsize=15)

    # Actual plotting
    ax.scatter(oi_halpha, oiii_hbeta, s=9, color='k', zorder=10)

    # read in Mappings III models and overplot
    mappings_oi_halpha_v100, mappings_oi_halpha_v125, mappings_oi_halpha_v150,\
    mappings_oi_halpha_v175, mappings_oi_halpha_v200, mappings_oi_halpha_v225,\
    mappings_oi_halpha_v250, mappings_oi_halpha_v300, mappings_oi_halpha_v500,\
    mappings_oi_halpha_v800,\
    mappings_nii_halpha_v100, mappings_nii_halpha_v125, mappings_nii_halpha_v150,\
    mappings_nii_halpha_v175, mappings_nii_halpha_v200, mappings_nii_halpha_v225,\
    mappings_nii_halpha_v250, mappings_nii_halpha_v300, mappings_nii_halpha_v350,\
    mappings_nii_halpha_v400, mappings_nii_halpha_v450, mappings_nii_halpha_v500,\
    mappings_nii_halpha_v800,\
    mappings_oiii_hbeta_v100, mappings_oiii_hbeta_v125, mappings_oiii_hbeta_v150,\
    mappings_oiii_hbeta_v175, mappings_oiii_hbeta_v200, mappings_oiii_hbeta_v225,\
    mappings_oiii_hbeta_v250, mappings_oiii_hbeta_v300, mappings_oiii_hbeta_v350,\
    mappings_oiii_hbeta_v400, mappings_oiii_hbeta_v450, mappings_oiii_hbeta_v500,\
    mappings_oiii_hbeta_v800,\
    mappings_sii_halpha_v100, mappings_sii_halpha_v125, mappings_sii_halpha_v150,\
    mappings_sii_halpha_v175, mappings_sii_halpha_v200, mappings_sii_halpha_v225,\
    mappings_sii_halpha_v250, mappings_sii_halpha_v300, mappings_sii_halpha_v500,\
    mappings_sii_halpha_v800 = bpt.mappings_oi_nii_sii()
    
    # Classification lines
    y_agn_hii_line = 1.33 + 0.73 / (np.arange(-2.5, -0.8, 0.01) + 0.59)
    y_liner_seyfert_line = 1.30 + 1.18 * np.arange(-1.12, 0, 0.01)

    ax.plot(np.arange(-2.5, -0.8, 0.01), y_agn_hii_line, '-', color='k')
    ax.plot(np.arange(-1.12, 0, 0.01), y_liner_seyfert_line, '--', color='k')

    ax.plot(mappings_oi_halpha_v125, mappings_oiii_hbeta_v125, '.-', lw=2, label='125 km/s')
    ax.plot(mappings_oi_halpha_v175, mappings_oiii_hbeta_v175, '.-', lw=2, label='175 km/s')
    ax.plot(mappings_oi_halpha_v200, mappings_oiii_hbeta_v200, '.-', lw=2, label='200 km/s')
    ax.plot(mappings_oi_halpha_v300, mappings_oiii_hbeta_v300, '.-', lw=2, label='300 km/s')
    ax.plot(mappings_oi_halpha_v500, mappings_oiii_hbeta_v500, '.-', lw=2, label='500 km/s')
    ax.plot(mappings_oi_halpha_v800, mappings_oiii_hbeta_v800, '.-', lw=2, label='800 km/s')

    ax.set_xlim(-2.0,0)
    ax.set_ylim(-1,1)

    # labels
    seyfertbox = TextArea('Seyfert', textprops=dict(color='k', size=16))
    anc_seyfertbox = AnchoredOffsetbox(loc=2, child=seyfertbox, pad=0.0, frameon=False,\
                                         bbox_to_anchor=(0.35, 0.93),\
                                         bbox_transform=ax.transAxes, borderpad=0.0)
    ax.add_artist(anc_seyfertbox) 

    linerbox = TextArea('LINER', textprops=dict(color='k', size=16))
    anc_linerbox = AnchoredOffsetbox(loc=2, child=linerbox, pad=0.0, frameon=False,\
                                         bbox_to_anchor=(0.8, 0.45),\
                                         bbox_transform=ax.transAxes, borderpad=0.0)
    ax.add_artist(anc_linerbox) 

    hiibox = TextArea('HII', textprops=dict(color='k', size=16))
    anc_hiibox = AnchoredOffsetbox(loc=2, child=hiibox, pad=0.0, frameon=False,\
                                         bbox_to_anchor=(0.2, 0.2),\
                                         bbox_transform=ax.transAxes, borderpad=0.0)
    ax.add_artist(anc_hiibox)

    # Add labels for velocities
    ax.legend(loc=6)

    # turn on minorticks
    ax.minorticks_on()

    # Save
    fig.savefig(taffydir + 'oi_halpha_bpt_xliners.eps', dpi=300, bbox_inches='tight')

    sys.exit(0)