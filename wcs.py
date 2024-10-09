"""
some tools for world coordinate systems

"""

from PyGalKin import *

import astropy.wcs as pywcs

def new_wcs(ra,dec,scale,pix=[1,1]):
    """
    ra, dec as float in dregrees
    scale in ""/pix
    assumes the image is already north-up east-left
    """
    scale /= 3600.
    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.equinox = 2000.0

    wcs.wcs.crpix = [ra, dec]
    wcs.wcs.cdelt = [scale, scale]
    wcs.wcs.crval = [0, -90]
    wcs.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    return wcs

def newer_wcs(ra,dec,scale,pix=[1,1]):
    """
    ra, dec as float in dregrees
    scale in ""/pix
    assumes the image is already north-up east-left
    """
    scale /= 3600.
    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.equinox = 2000.0

    wcs.wcs.crval = [ra, dec]
    wcs.wcs.cdelt = [scale, scale]
    wcs.wcs.crpix= pix
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    return wcs
