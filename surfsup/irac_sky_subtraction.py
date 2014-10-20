#!/usr/bin/env python

import numpy as np
from scipy import ndimage, signal
import pyfits
import os

"""
Try different tricks to subtract intra-cluster light (ICL).
"""

class IRACSkySubtrac(object):
   def __init__(self, imagename):
      self.imagename = imagename
      self.image = pyfits.getdata(imagename)
      self.hdr = pyfits.getheader(imagename)
      self.skyimage = np.zeros(self.image.shape)

   def median_filter(self, skyimagename, newimagename, size=51):
      """
      Performs a median filtering to estimate ICL, and then subtract from the 
      original image.
      """
      self.skyimage = signal.medfilt2d(self.image, (size, size))
      if os.path.exists(skyimagename):
         os.remove(skyimagename)
      if os.path.exists(newimagename):
         os.remove(newimagename)
      pyfits.append(skyimagename, self.skyimage, self.hdr)
      pyfits.append(newimagename, self.image-self.skyimage, self.hdr)
