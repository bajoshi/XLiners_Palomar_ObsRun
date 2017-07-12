from astropy.io import fits
import sys

def fix_la_cosmic_header():

    ext_dir = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/'
    date = '20170509/'
    filename = 'tspec0062'
    
    orig_name = ext_dir + date + filename + '.fits'
    clean_name = ext_dir + date + filename + '_cleaned.fits'
    
    orig = fits.open(orig_name)
    cleaned = fits.open(clean_name)
    
    cleaned[0].header = orig[0].header
    cleaned.writeto(clean_name, clobber=True)  # overwrites cleaned file with proper header
    
    return None

def replace_header_for_setjd():

    # MAKE SURE that you've run setjd on the flat fielded image before running this funciton

    ext_dir = '/Volumes/Bhavins_backup/ipac/Palomar_data/2017/'
    date = '20170509/'
    filename = 'tspec0062'
    
    flat_name = ext_dir + date + filename + '_f.fits'
    clean_name = ext_dir + date + filename + '_cleaned.fits'

    flat = fits.open(flat_name)
    cleaned = fits.open(clean_name)

    cleaned[0].header = flat[0].header
    cleaned.writeto(clean_name, clobber=True)  # overwrites cleaned file with proper header

    return None

if __name__ == '__main__':

    # function names tell you what the functions do. Uncomment the line you want.
    # The two functions will be combined at some point

    #fix_la_cosmic_header()
    replace_header_for_setjd()
    
    sys.exit(0)