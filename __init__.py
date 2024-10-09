"""
Doc-String for the __init__.py file.
look into the doc string for the PyCigale submodule
for more information.
"""

__author__ = 'T. Marquart, J. Staabis, A. Zien'
__version__= 0.2

__PAR__= {}

# scipy etc.
# (pylab is imported in plot.py)
import numpy as N
import astropy.io.fits as pyfits
import scipy as S
from scipy import interpolate
from scipy import signal as Sig
from scipy.optimize import leastsq as LS
import numpy.ma as ma
from numpy.ma import masked_where,masked_array,mask_or
from decimal import Decimal
from scipy.fftpack import fft
#from filtfilt import filtfilt
from mpfit import mpfit
from math import pi,e,radians,degrees
import matplotlib.mlab as mlab
import astropy.wcs as pywcs

# OS etc.
import os
import os.path as path
import sys
from time import sleep
from copy import copy
#import wx
import pickle
from sqlite3 import dbapi2 as sqlite
#from sdss import sqlcl

##############
# own things #
##############
from .tool import *
from . import wcs as WCS
from . import PyCigale as C
#import PyArgus as A
from . import InOutput as IO
#import ModelVF as MV
#import Gauss as G
from . import fit as F
from . import plot as P
#import MacroGen as MG
from . import db as DB
from . import interact as I

# Good ol' IDL
#try: from pyIDL import idl as IDL
#except: print 'Could not import IDL'



