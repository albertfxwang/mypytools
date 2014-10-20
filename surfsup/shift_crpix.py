#!/usr/bin/env python

import numpy as np
import pyfits
import os, sys
import pywcs

# TFIT requires that CRPIX's are half-integers, so here we go...
# I will subtract CRPIX's of each image until they are half-integers, 
# and record the subtracted values in the image header
# This will also require an update of the CRVALs to be consistent with the 
# previous WCS

img_ext = ['drz','cov','std','unc']

def shift_crpix(drzimage):
   root = os.path.splitext(drzimage)[0][:-3]
   for ext in img_ext:
      fname = root + ext + '.fits'
      if os.path.exists(fname):
         h = pyfits.open(fname, mode='update')
         crpix1 = h[0].header['crpix1']
         crpix2 = h[0].header['crpix2']
         q1 = int(crpix1 / 0.5)
         q2 = int(crpix2 / 0.5)
         # q1, q2 equal to the nearest smaller half-integer divided by 0.5
         # so q1 and q2 should be odd numbers; if not, subtract by 1
         if q1 % 2 == 0:
            q1 = q1 - 1
         if q2 % 2 == 0:
            q2 = q2 - 1
         # t1, t2 are the closest smaller half-integers
         t1 = q1 * 0.5
         t2 = q2 * 0.5
         resid1 = crpix1 - t1
         resid2 = crpix2 - t2
         # Now re-calculate CRVALs
         wcs_img = pywcs.WCS(h[0].header)
         crval1, crval2 = wcs_img.wcs_pix2sky([[t1, t2]], 1)[0]
         h[0].header['crpix1'] = t1
         h[0].header['dcrpix1'] = resid1
         h[0].header['crpix2'] = t2
         h[0].header['dcrpix2'] = resid2
         crval1_0 = h[0].header['crval1']
         crval2_0 = h[0].header['crval2']
         h[0].header['crval1_0'] = crval1_0
         h[0].header['crval2_0'] = crval2_0
         h[0].header['crval1'] = crval1
         h[0].header['crval2'] = crval2
         h.flush()
         h.close()
         print "%s processed." % fname

