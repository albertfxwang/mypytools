#!/usr/bin/env python

# Match the image footprints... started from segmap_30mas.py

import numpy as np
import os, sys
import pywcs, pyfits

def match_images(inputimage, refimage, output=None):
# def match_images(segimage, drzimage_v2, output=None):
   """
   Match the footprint of inputimage to that of refimage.
   """
   hdr1 = pyfits.getheader(inputimage)
   hdr2 = pyfits.getheader(refimage)
   wcs1 = pywcs.WCS(hdr1)
   wcs2 = pywcs.WCS(hdr2)
   sky00 = wcs2.wcs_pix2sky([[1, 1]], 1)
   corner00 = wcs1.wcs_sky2pix(sky00, 1)[0]
   corner00 = np.around(corner00).astype('int')
   nx2 = hdr2['naxis1']
   ny2 = hdr2['naxis2']
   sky11 = wcs2.wcs_pix2sky([[nx2, ny2]], 1)
   corner11 = wcs1.wcs_sky2pix(sky11, 1)[0]
   corner11 = np.around(corner11).astype('int')
   xlo1, ylo1 = corner00
   xhi1, yhi1 = corner11
   print "[xlo:xhi, ylo:yhi]", xlo1, xhi1, ylo1, yhi1
   if output == None:
      return xlo1, xhi1, ylo1, yhi1
   newshape = (xhi1-xlo1+1, yhi1-ylo1+1)
   assert newshape == (nx2, ny2), "Shape of new seg map does not match the shape of the input drz image..."
   # Now make a cutout of the seg map
   if os.path.exists(output):
      os.remove(output)
   im1 = pyfits.getdata(inputimage)
   im1_new = seg[ylo1:yhi1+1, xlo1:xhi1+1]
   # hdr1['crpix1'] = hdr1['crpix1'] - xlo1
   # hdr1['crpix2'] = hdr1['crpix2'] - ylo1
   hdr1['crpix1'] = hdr2['crpix1']
   hdr1['crpix2'] = hdr2['crpix2']
   print "shape of GOODS/ACS v2 mosaic: [%d, %d]" % (hdr2['naxis1'], hdr2['naxis2'])
   print "new shape of the segmap: [%d, %d]" % (im1_new.shape[1], im1_new.shape[0])
   pyfits.append(output, im1_new, hdr1)