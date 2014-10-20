#!/usr/bin/env python

import numpy as np
import pyfits
import os, sys

"""
Some IRAF/PyRAF tasks rewritten to use only python and no IRAF.
Pixel indices start at 1, not 0, to comply with IRAF conventions.
"""

def imcopy(input_file, section, output_file):
   if os.path.exists(output_file):
      raise IOError,  "%s already exists." % output_file
   try: 
      xmin, xmax, ymin, ymax = section
   except ValueError:
      raise ValueError, "Please give a list or array for the second argument."
   hdr = pyfits.getheader(input_file)
   hdr['crpix1'] = hdr['crpix1'] - (xmin-1)
   hdr['crpix2'] = hdr['crpix2'] - (ymin-1)
   input_image = pyfits.getdata(input_file)
   xmin = np.maximum(xmin, 1)
   ymin = np.maximum(ymin, 1)
   print xmin, xmax, ymin, ymax
   new_image = input_image[ymin-1:ymax,xmin-1:xmax]
   pyfits.append(output_file, new_image, hdr)
