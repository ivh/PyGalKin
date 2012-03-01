"""
some tools for world coordinate systems

"""

from PyGalKin import *

import pywcs

def new_wcs(ra,dec,scale):
    """
    ra, dec as float in dregrees
    scale in ""/pix
    assumes the image is already north-up east-left
    """
    scale /= 3599.
    wcs = pywcs.WCS(naxis=2)
    wcs.wcs.crpix = [ra, dec]
    wcs.wcs.cdelt = [scale, scale]
    wcs.wcs.crval = [0, -90]
    wcs.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    wcs.wcs.set_pv([(2, 1, 45.0)])
    wcs.wcs.equinox = 2000.0
    return wcs
