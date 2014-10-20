#!/usr/bin/env python

from numpy import *
import pyfits

# zero-points for the GOODS/ACS v2 mosaics
zpt_goodsv2 = {'f435w':25.65288, 
               'f606w':26.49341, 
               'f775w':25.64053, 
               'f850lp':24.84315}

# zero-points for the HUDF/ACS mosaics
zpt_udf = {'f435w':25.673, 
           'f606w':26.486, 
           'f775w':25.654, 
           'f850lp':24.862}

# zero-points for the CANDELS mosaics (at least in GOODS-S)
zpt_candels = {'vimos_u':26.158,
               'acs_f435w':25.673,
               'acs_f606w':26.505, 
               'acs_f775w':25.678,
               'acs_f814w':25.947,
               'acs_f850lp':24.867,
               'wfc3_f098m':25.667,
               'wfc3_f105w':26.270,
               'wfc3_f125w':26.250,
               'wfc3_f160w':25.960}


# HST zeropoint-related calculations

def wfc3_abmag_zeropoint(photflam, photplam):
   # calculate zeropoint in ABmag using PHOTFLAM and PHOTPLAM keywords in WFC3 images
   abmag = -2.5 * log10(photflam) - 21.10 - 5.0 * log10(photplam) + 18.6821
   return abmag


def wfc3_stmag_zeropoint(photflam):
   # calculate zeropoint in STmag using PHOTFLAM and PHOTZPT keywords in WFC3 images
   stmag = -2.5 * log10(photflam) - 21.10
   return stmag

def get_abmag_zpt(image):
  hdr = pyfits.getheader(image)
  photflam = hdr['photflam']
  photplam = hdr['photplam']
  return wfc3_abmag_zeropoint(photflam, photplam)

