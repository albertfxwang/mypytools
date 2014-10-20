#!/usr/bin/env python

import numpy as np
from pyraf import iraf
import pyfits
import pywcs
import os


def cutout(image, ra, dec, xwidth, ywidth, scale):
   """
   Make a cutout around position (ra, dec) with the dimensions of the cutout
   being (xwidth, ywidth).
   image: file name of the image
   ra, dec: sky coordinates of the center of cutout, in degrees
   xwidth, ywidth: dimensions of the cutout, **in arcsec**
   scale: pixel scale of the input image, **in arcsec**
   """
   hdr = pyfits.getheader(image)
   # First calculate the image coordinates, based on 1-indexing
   wcs = pywcs.WCS(hdr)
   xypos = wcs.wcs_sky2pix([[ra, dec]], 1)[0]
   xc = int(round(xypos[0]))
   yc = int(round(xypos[1]))
   xw_pix = int(round(xwidth / (2. * scale)))
   yw_pix = int(round(ywidth / (2. * scale)))
   xmin = xc - xw_pix
   xmax = xc + xw_pix
   ymin = yc - yw_pix
   ymax = yc + yw_pix
   outname = raw_input('Please give the file name of the output image: ')
   if not len(outname):
      outname = 'cutout.fits'
      print "No file name is given; use cutout.fits"
   if os.path.exists(outname):
      os.remove(outname)
   iraf.imcopy(image+'[%d:%d,%d:%d]' % (xmin, xmax, ymin, ymax), outname)