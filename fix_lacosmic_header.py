from astropy.io import fits
import sys

def fix_la_cosmic_header(wd, filename, for_setjd=False):
    
    if for_setjd:
        orig_name = wd + filename + '_f.fits'
        clean_name = wd + filename + '_cleaned.fits'
    else:
        orig_name = wd + filename + '.fits'
        clean_name = wd + filename + '_cleaned.fits'
    
    orig = fits.open(orig_name)
    cleaned = fits.open(clean_name)
    
    cleaned[0].header = orig[0].header
    cleaned.writeto(clean_name, clobber=True)  # overwrites cleaned file with proper header
    
    return None

if __name__ == '__main__':

    # function names tell you what the functions do. Uncomment the line you want.
    # The two functions will be combined at some point

    fix_la_cosmic_header()
    
    sys.exit(0)