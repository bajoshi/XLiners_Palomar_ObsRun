from astropy.io import fits

orig_name = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/20170511/tspec0026.fits'
clean_name = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/20170511/tspec0026_cleaned.fits'

orig = fits.open(orig_name)
cleaned = fits.open(clean_name)

cleaned[0].header = orig[0].header
cleaned.writeto(clean_name, clobber=True)  # overwrites cleaned file with proper header