#!/usr/bin/env python

import numpy as np
import pyfits, pywcs
from pyraf import iraf
import os, sys


class cutouts(object):
   """
   A base class for making cutouts of the same object in multiple filters.
   """
   def __init__(self, images, filters):
      """
      filters is a list of strings representing filter names.
      """
      self.images = {}
      self.filters = filters
      self.wcs = {}
      self.pixscale = {}  # in arcsec
      for i in range(len(filters)):
         print filters[i]
         self.images[filters[i]] = images[i]
         h = pyfits.getheader(images[i])
         if 'CPDIS1' in h.keys():
            print "Keyword CPDIS1 exists in header, and it might create errors for pywcs..."
         try:
            self.wcs[filters[i]] = pywcs.WCS(h, naxis=2)
         except:
            print "WCS for %s is not supported..." % filters[i].upper()
         scale = np.sqrt(h['cd1_1']**2 + h['cd1_2']**2) * 3600.
         scale = np.round(scale, 3)
         print "Pixel scale in %s is %.3f arcsec" % (filters[i], scale)
         self.pixscale[filters[i]] = scale

   def cut_xy(self, x, y, band, width, name, norm=False):
      # width is in arcsec
      imgname = '%s_%s.fits' % (name, band)
      if os.path.exists(imgname):
         os.remove(imgname)
      width_pix = width / self.pixscale[band]
      xmin = int((x - width_pix / 2.))
      xmax = int((x + width_pix / 2.))
      ymin = int((y - width_pix / 2.))
      ymax = int((y + width_pix / 2.))
      # Enforce that the image have equal number of pixels in both dimensions
      npix = np.min([xmax-xmin, ymax-ymin])
      if xmax - xmin > npix:
         xmin += 1
      elif ymax - ymin > npix:
         ymin += 1
      self.xmin = xmin
      self.ymin = ymin
      iraf.imcopy('%s[%d:%d,%d:%d]' % (self.images[band], xmin, xmax, ymin, ymax),
                  imgname, verbose=False)
      if norm:
         imgsum = pyfits.getdata(imgname).ravel().sum()
         imgsum = np.abs(imgsum)
         iraf.imarith(imgname, '/', imgsum, 'temp.fits')
         os.remove(imgname)
         os.system('mv temp.fits %s' % imgname)
      return imgname

   def cut_xy_all(self, x, y, bands, width, name, norm=False):
      imgnames = map(lambda b: self.cut_xy(x, y, b, width, name, norm=norm), 
                     self.filters)
      self.display(imgnames)

   def cut_radec(self, ra, dec, band, width, name, norm=False):
      # width is in arcsec
      x, y = self.wcs[band].wcs_sky2pix(np.array([[ra, dec]]), 1)[0]
      self.x_image = x
      self.y_image = y
      return self.cut_xy(x, y, band, width, name, norm=norm)

   def cut_radec_all(self, ra, dec, bands, width, name, norm=False):
      imgnames = map(lambda b: self.cut_radec(ra, dec, b, width, name, norm=norm), 
                     self.filters)
      self.display(imgnames)

   def display(self, imgnames):
      cmd = 'ds9 -multiframe '
      for img in imgnames:
         cmd += img + ' -zoom to fit -scale zscale '
      cmd += '-match frame wcs & '
      print cmd
      os.system(cmd)

class Cutouts(cutouts):
   # A name proxy for cutout class (so I don't need to rename it)
   def __init__(self, images, filters):
      super(Cutouts, self).__init__(images, filters)