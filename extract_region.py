import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS

if __name__ == '__main__':

    #def sedigism_extract():
        # Inputs of function
        line = 13   # V0: 13/18, V1: 13/18/'13'/'18'/'13CO'/'C18O' all supported. note the filename in the webpage is lowercase as c18o21
        # fieldcenter = '001'   # should be calculated automatically from the input.
        clon = 355.7  # Galatic longtitude of extracted cube [deg]
        clat = -0.4 # Galatic latitutde of extracted cube [deg]
        extSizeL = 900 # half size of the extracted cube [arcsecond] (e.g. +/- 100")
         

        # Read in datacube
        fitfile = 'G'+fieldcenter+'_'+line+'_Tmb_DR1.fits'
        datacube = fits.open(fitfile)[0]

        # Convert subcube center from world to pixel coordinate
        cworld = SkyCoord(l=clon,b=clat,unit=(u.deg, u.deg),frame='galactic')
        w = WCS(datacube.header)
        convert = w.world_to_pixel(cworld,3000 * u.m / u.s)
        cpix_long = convert[0]
        cpix_lat  = convert[1]


        # Convert subcube size from world to pixel coordinate
        pix = 9.5 # pixel size [arcseconds]
        extsize_pix = extsize/pix

        # Determine the range of datacube to extract
        long_range = np.round(np.array([cpix_long-extsize_pix, cpix_long+extsize_pix]),decimals=0).astype(int)
        lat_range  = np.round(np.array([cpix_lat-extsize_pix, cpix_lat+extsize_pix]),decimals=0).astype(int)

        # Extract subcube and save to fitsfile
        subcube = datacube.data[:,lat_range[0]:lat_range[1]+1,long_range[0]:long_range[1]+1]
        subcube_header = datacube.header.copy()
        subcube_header['NAXIS1']  = np.shape(subcube)[2]
        subcube_header['NAXIS2']  = np.shape(subcube)[1]
        subcube_header['DATAMIN'] = np.nanmin(subcube)
        subcube_header['DATAMAX'] = np.nanmax(subcube)
        subcube_header['CRVAL1']  = clon
        subcube_header['CRPIX1']  = cpix_long - long_range[0]
        subcube_header['CRVAL2']  = clat
        subcube_header['CRPIX2']  = cpix_lat - lat_range[0]

        fits.writeto('long_lat_size_Tmb.fits',subcube,header=subcube_header,overwrite=True)
