#!/usr/bin/env python

import numpy as np
import os
from ImageTools import images
from pygoods import sextractor

class InputXY(object):
   """
   Calculate input (x, y) positions.
   """
   def __init__(self, drzimage, flagimage):
      # define input science and flag images
      self.drzimage = images.FITSImage(drzimage)
      self.flgimage = images.FITSImage(flagimage)
      self.xinput = []
      self.yinput = []

   def reset(self):
      self.xinput = []
      self.yinput = []

   def inputxy_mindist(self, ngals, catalog, xcol='x_image', 
                       ycol='y_image', mindist=0.4, pixscale=0.06,
                       maxflag=1):
      """
      Calculate ngals input positions at least mindist (in arcsec) away from 
      all detected sources in catalog. Pixel scale is pixscale (in arcsec).
      Write the generate (x, y) position to output.
      """
      c = sextractor(catalog)
      xsrc = getattr(c, xcol)
      ysrc = getattr(c, ycol)
      nlines, ncols = self.drzimage.data.shape
      n = 0
      mindist_pix = mindist / pixscale
      while n < ngals:
         x = np.random.uniform(0, ncols)
         y = np.random.uniform(0, nlines)
         xr = int(np.floor(x))
         yr = int(np.floor(y))
         mindist_n = np.sqrt((x-xsrc)**2 + (y-ysrc)**2).min()
         if (mindist_n >= mindist_pix) & (self.flgimage.data[yr,xr]<maxflag):
            self.xinput += [x]
            self.yinput += [y]
            n += 1
      self.xinput = np.array(self.xinput)
      self.yinput = np.array(self.yinput)
      return self.xinput, self.yinput

   def write_xy(self, output):
      f = open(output, 'wb')
      f.write('# 1 X\n')
      f.write('# 2 Y\n')
      for i in range(len(self.xinput)):
         f.write('%.2f  %.2f  ' % (self.xinput[i], self.yinput[i]))
         f.write('\n')
      f.close()

   def write_reg(self, filename='inputxy.reg'):
      f = open(filename, 'wb')
      for i in range(len(self.xinput)):
         f.write('image; circle(%.2f, %.2f, 1")\n' % (self.xinput[i], self.yinput[i]))
      f.close()

