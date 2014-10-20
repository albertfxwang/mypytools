#!/usr/bin/env python

import numpy as np
import pyfits
from pygoods import sextractor, Ftable
import os
from pyraf import iraf
iraf.stsdas()

def make_covar_map(idmap, catalog):
   """
   Given the IDMAP image, replace the ID by the covariance index of the object,
   in order to visualize the distribution of covariance indices in the image.
   """
   c = sextractor(catalog)
   cvdict = {}
   for i in range(len(c)):
      cvdict[c.objectid[i]] = c.maxcvratio[i]

   # Now make covariance id map
   covar_idmap = 'covar_' + idmap
   if os.path.exists(covar_idmap):
      os.remove(covar_idmap)
   iraf.imcalc(idmap, covar_idmap, 'im1*1.0', pixtype='double')
   #os.system('cp %s %s' % (idmap, covar_idmap))
   h = pyfits.open(covar_idmap, mode='update')
   data = h[0].data
   datacopy = np.copy(data)
   for k, v in cvdict.iteritems():
      datacopy[data==k] = v
   h[0].data = datacopy
   h.flush()
   h.close()