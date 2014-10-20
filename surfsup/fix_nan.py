#!/usr/bin/env python

import numpy as np
import pyfits

def fix_nan(image, replace=0.):
   """
   Fix the NAN's by replacing them with a constant.
   Usually, replace NAN by 0 if it's a science exposure or a weight image, 
   and replace NAN by a huge number (like 1.e10) if it's an RMS image.
   """
   h = pyfits.open(image, mode='update')
   imgdata = h[0].data
   imgdata = np.where(np.isnan(imgdata), replace, imgdata)
   h[0].data = imgdata
   h.flush()
   h.close()

   
