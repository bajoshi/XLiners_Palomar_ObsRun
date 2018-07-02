from __future__ import division

from astropy.io import fits
import numpy as np

import os
import sys

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

home = os.getenv('HOME')  # Does not have a trailing slash at the end
xliners_dir = home + '/Desktop/ipac/XLiners_Palomar_ObsRun/'

sys.path.append(xliners_dir)

if __name__ == '__main__':

    # define directories for both nights
    wd_n1 = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/fluxcal_workdir2018_2017Anight1/'
    wd_n2 = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/fluxcal_workdir2018_2017Anight3/'
    
    # read all fits files needed
    n1a = fits.open(wd_n1 + 'hip64248_AB_bksub.fits')
    n1b = fits.open(wd_n1 + 'hip64248_BA_bksub.fits')

    n2a = fits.open(wd_n2 + 'hip64248_AB_bksub.fits')
    n2b = fits.open(wd_n2 + 'hip64248_BA_bksub.fits')

    # Assign arrays and normalize
    # The normalization will not affect the shape of hte spectrum.
    # It will only change the scale. And we only want to compare
    # the shape in any case. This is because the scale is likely
    # to be different for each night given that the light received,
    # in terms of absolute counts, even though from the same object 
    # might be different due to observing conditions like seeing
    # and airmass.
    # All A positions
    ja_n1 = n1a[0].data[0][0] #/ np.nanmedian(n1a[0].data[0][0])
    ha_n1 = n1a[0].data[0][1] #/ np.nanmedian(n1a[0].data[0][1])
    ka_n1 = n1a[0].data[0][2] #/ np.nanmedian(n1a[0].data[0][2])
    ja_n2 = n2a[0].data[0][0] #/ np.nanmedian(n2a[0].data[0][0])
    ha_n2 = n2a[0].data[0][1] #/ np.nanmedian(n2a[0].data[0][1])
    ka_n2 = n2a[0].data[0][2] #/ np.nanmedian(n2a[0].data[0][2])

    # All B positions
    jb_n1 = n1b[0].data[0][0] #/ np.nanmedian(n1b[0].data[0][0])
    hb_n1 = n1b[0].data[0][1] #/ np.nanmedian(n1b[0].data[0][1])
    kb_n1 = n1b[0].data[0][2] #/ np.nanmedian(n1b[0].data[0][2])
    jb_n2 = n2b[0].data[0][0] #/ np.nanmedian(n2b[0].data[0][0])
    hb_n2 = n2b[0].data[0][1] #/ np.nanmedian(n2b[0].data[0][1])
    kb_n2 = n2b[0].data[0][2] #/ np.nanmedian(n2b[0].data[0][2])

    # Plot to compare
    gs = gridspec.GridSpec(10,21)
    gs.update(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.0, hspace=0.0)

    fig = plt.figure()

    # All As; top panels
    ax1 = fig.add_subplot(gs[:5,:7])
    ax2 = fig.add_subplot(gs[:5,7:14])
    ax3 = fig.add_subplot(gs[:5,14:])

    # All Bs; bottom panels
    ax4 = fig.add_subplot(gs[5:,:7])
    ax5 = fig.add_subplot(gs[5:,7:14])
    ax6 = fig.add_subplot(gs[5:,14:])

    # Jband A pos
    ax1.plot(np.arange(len(ja_n1)), ja_n1)
    ax1.plot(np.arange(len(ja_n2)), ja_n2)

    # Hband A pos
    ax2.plot(np.arange(len(ha_n1)), ha_n1)
    ax2.plot(np.arange(len(ha_n2)), ha_n2)

    # Kband A pos
    ax3.plot(np.arange(len(ka_n1)), ka_n1)
    ax3.plot(np.arange(len(ka_n2)), ka_n2)

    # Jband B pos
    ax4.plot(np.arange(len(jb_n1)), jb_n1)
    ax4.plot(np.arange(len(jb_n2)), jb_n2)

    # Hband B pos
    ax5.plot(np.arange(len(hb_n1)), hb_n1)
    ax5.plot(np.arange(len(hb_n2)), hb_n2)

    # Kband B pos
    ax6.plot(np.arange(len(kb_n1)), kb_n1)
    ax6.plot(np.arange(len(kb_n2)), kb_n2)

    plt.show()

    plt.cla()
    plt.clf()
    plt.close()

    # Divide 
    ja_res = ja_n1 / ja_n2
    ha_res = ha_n1 / ha_n2
    ka_res = ka_n1 / ka_n2

    fig2 = plt.figure()
    ax1 = fig2.add_subplot(gs[:,:7])
    ax2 = fig2.add_subplot(gs[:,7:14])
    ax3 = fig2.add_subplot(gs[:,14:])

    ax1.plot(np.arange(len(ja_res)), ja_res)
    ax2.plot(np.arange(len(ha_res)), ha_res)
    ax3.plot(np.arange(len(ka_res)), ka_res)

    plt.show()

    # Close fits HDUs and exit
    n1a.close()
    n1b.close()
    n2a.close()
    n2b.close()

    sys.exit(0)