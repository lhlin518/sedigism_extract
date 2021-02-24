import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS

if __name__ == '__main__':

    #def sedigism_extract( coord, size, line=13, inputPath=None, outputPath=None, outputFilename=None ):
        # Inputs of function
        # In astronomy, lonLeft > lonRight, latBottom < latTop
        lonLeft = 355.7     # clon + extSizeL / 2
        lonRight = 354.9    # clon - extSizeL / 2
        latBottom = -0.2    # clat - extSizeB / 2
        latTop  = +0.1      # clat + extSizeB / 2
        clon = ( lonLeft + lonRight ) / 2
        clat = ( latBottom + latTop ) / 2
        extSizeL = ( lonLeft - lonRight ) / 2   # in deg, or we use quantities here.
        extSizeB = ( latBottom - latTop ) / 2
        
        
        line = 13   
        
        if line == 13 or line == '13' or line == '13CO':
            lineFilename = '13CO21'     # the filename in the webpage
            lineHeader = '13CO J=2-1'   # the name in header['COMMENT'] or somewhere
            freqLine = 220398.6842 * u.MHz
        elif line == 18 or line == '18' or line == 'C18O':
            lineFilename = 'C18O21'     # the filename in the webpage
            lineHeader = 'C18O J=2-1'   # the name in header['COMMENT'] or somewhere
            freqLine = 219560.3541 * u.MHz
        else:
            sys.exit("Please input line as \'13CO\' or \'C18O\'." )
        
        # fieldcenter = '001'   # should be calculated automatically from the input.
        #clon = 355.7  # Galatic longtitude of extracted cube [deg]
        #clat = -0.4 # Galatic latitutde of extracted cube [deg]
        #extSizeL = 900 # half size of the extracted cube [arcsecond] (e.g. +/- 100")
        
        # if the longitude range covers two integer, it is too large.
        if np.floor( lonLeft ) - np.ceil( lonRight ) > 0:
            sys.exit('This script cannot extract data from two blocks.')
        
        lonBlockCenter = np.round(clon).astype(int)
        filename = 'G{:d}_'.format(lonBlockCenter) + lineFilename + '_Tmb_DR1.fits'
        
        # search the path
        if inputPath == None:
            inputPath = './'         # default to current folder
        elif inputPath == 'webpage' or inputPath == 'Web':
            inputPath = 'web'
        elif inputPath[-1] != '/':
            inputPath = inputPath + '/'   # end with '/'
        weblink = 'https://sedigism.mpifr-bonn.mpg.de/SEDIGISM_DATACUBES_DR1c/'
        for path in [ inputPath, './', 'MPIfR', 'web' ]:
            if path != 'MPIfR' and path != 'web':
                try:
                    f = open( path + filename )
                    f.close()
                    break
                except IOError:
                    print( path + ' does not cotain ' + filename )
                    continue
            elif path == 'MPIfR':
                try:
                    path = '/aux/atlasgalb/www-sedigism/SEDIGISM_DATACUBES_DR1c/'
                    f = open( path + filename )
                    f.close()
                    break
                except IOError:
                    print( path + ' does not cotain ' + filename )
                    continue
            elif path == 'web':
                path == './'
                with requests.get( weblink + filename, timeout=5, stream=True ) as req:
                    with open( './' + filename, 'wb' ) as f:   
                        for chunk in req.iter_content(chunk_size=1024*1024):
                            if chunk:
                                f.write( chunk )
                break
        # End getdata
        # open the fits
        data, header = fits.getdata( path + filename, header=True )

        
        
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
