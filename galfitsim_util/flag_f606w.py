#!/usr/bin/env python

import numpy as np
from pygoods import *
import os, sys
import pyfits

"""
To make F606W flag image in UDF, to mask out the detections in F606W that are
not picked up in F160W, to avoid those causing trouble for GALFIT.
"""

def flag_f606w(flag0, seg_vdetect, cat_hdetect, cat_vdetect, tol=1.0):
   # tol: tolerance in arcsec for matching sources between the H-detected and
   # V-detected catalogs
   ch = sextractor(cat_hdetect)
   cv = sextractor(cat_vdetect)
   seg = pyfits.getdata(seg_vdetect)
   flag = pyfits.getdata(flag0)
   flagged_id = []
   for i in range(len(cv)):
      angdist = angsep.angsep(cv.alpha_j2000[i], cv.delta_j2000[i],
                              ch.alpha_j2000, ch.delta_j2000)
      if angdist.min() * 3600. > tol:
         flagged_id += [cv.number[i]]

   flagged_id = np.array(flagged_id)
   print "len(flagged_id)", len(flagged_id)
   # Now make a new flag map for F606W
   flagnew = np.where((flag.ravel()>0) | (np.in1d(seg, flagged_id)), 1, 0)
   flagnew = flagnew.reshape(flag.shape).astype('int32')
   flagnewname = os.path.splitext(flag0)[0] + 'new.fits'
   if os.path.exists(flagnewname):
      os.remove(flagnewname)
   os.system('cp %s %s' % (flag0, flagnewname))
   h = pyfits.open(flagnewname, mode='update')
   h[0].data = flagnew
   h.flush()
   h.close()
