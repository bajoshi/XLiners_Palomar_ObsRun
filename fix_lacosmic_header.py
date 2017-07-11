from astropy.io import fits

ext_dir = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/'
date = '20170509/'
filename = 'xl53_BA_median'

orig_name = ext_dir + date + filename + '.fits'
clean_name = ext_dir + date + filename + '_cleaned.fits'

orig = fits.open(orig_name)
cleaned = fits.open(clean_name)

cleaned[0].header = orig[0].header
cleaned.writeto(clean_name, clobber=True)  # overwrites cleaned file with proper header