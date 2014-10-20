#!/usr/bin/env python

import numpy as np
import pyfits
import pywcs
from pyraf import iraf
from cutouts import cutouts
import os, sys
from pygoods import sextractor
import image_moments  # note that this uses the python 0-based indexing!!
from scipy.ndimage import interpolation, measurements


class PSFStack(cutouts):
   def __init__(self, catalog, image, band, width, norm=True, save=False):
      """
      The catalog (in SExtractor format) needs to have the following three columns:
      NUMBER
      RA
      DEC
      so that we can figure out where the sources are.
      filter: the name of the filter
      width: the width of the cutouts (in arcsec)
      """
      cutouts.__init__(self,  [image], [band])
      self.width = width
      self.c = sextractor(catalog)
      self.band = band
      self.norm = norm
      self.save = save
      print "There are a total of %d stars in the list." % len(self.c)

   def cut_all(self):
      # First clean all the obj*_cut.fits
      os.system('rm obj*.fits')
      self.xmin_list = []
      self.ymin_list = []
      self.x_image_list = []
      self.y_image_list = []
      for i in range(len(self.c)):
         self.cut_radec(self.c.ra[i], self.c.dec[i], self.band, self.width,
                        'obj%d' % (self.c.number[i]), norm=self.norm)
         self.xmin_list += [self.xmin]
         self.ymin_list += [self.ymin]
         self.x_image_list += [self.x_image]
         self.y_image_list += [self.y_image]

   def oversample_all(self, factor=10., interpolant='sinc'):
      """
      For each star, perform centroid finding and over-sample the images by 
      a factor supplied in the method call.
      """
      self.factor = factor
      self.interpolant = interpolant
      f = open('psflist', 'wb')
      for i in range(len(self.c)):
         imgname = 'obj%d_%s.fits' % (self.c.number[i], self.band)
         starimg = pyfits.getdata(imgname)
         pm = image_moments.moments(starimg)
         pm.firstorder()
         mag = 1./factor   # for use in iraf.imlintran
         # Now over-sample the image
         shape0 = starimg.shape
         shape1 = (np.array(shape0) * factor).astype('int')
         imgname_oversample = os.path.splitext(imgname)[0] + '_mag%.2f.fits' % (mag)
         # xin = self.x_image_list[i] - self.xmin_list[i]
         # yin = self.y_image_list[i] - self.ymin_list[i]
         xin = pm.y1 
         yin = pm.x1 
         iraf.imlintran(input=imgname,
                        output=imgname_oversample,
                        xrotation=0., yrotation=0.,
                        xmag=mag, ymag=mag, 
                        xin=xin, yin=yin,
                        # xin=pm.x1+1, yin=pm.y1+1,
                        # xout=shape1[0]/2., yout=shape1[1]/2.,
                        ncols=shape1[1], nlines=shape1[0],
                        interpolant=interpolant)
         f.write('%s\n' % imgname_oversample)
      f.close()

   def stack_all(self, mode='median', save=False):
      """
      Stack all PSF stars to produce a PSF image in the original pixel scale.
      """
      output = 'psf_stacked_raw_%s_%s_mag%.2f.fits' % (self.band, mode, 1./self.factor)
      if os.path.exists(output):
         os.remove(output)
      iraf.imcombine(input='@psflist', output=output, combine=mode,
                     reject='sigclip', lsigma=3.0, hsigma=3.0,
                     mclip='yes', nkeep=3)
      shape2 = np.array(pyfits.getdata(output).shape)
      # Now re-sample the image back to the original pixel scale
      output2 = 'psf_stacked_%s_%s.fits' % (self.band, mode)
      if os.path.exists(output2):
         os.remove(output2)
      # pm = image_moments.moments(pyfits.getdata(output))
      # pm.firstorder()
      pm1 = measurements.center_of_mass(pyfits.getdata(output))
      print "Oversampled PSF shape:", shape2
      print "Centroid before shift: %.2f, %.2f" % tuple(pm1)
      print "Shifted PSF image:", output2
      # AND REMEMBER TO SWAP X & Y if using IRAF!
      # iraf.imlintran(input=output, 
      #                output=output2,
      #                xrotation=0., yrotation=0.,
      #                xin=pm.y1, yin=pm.x1,
      #                xmag=self.factor, ymag=self.factor,
      #                ncols=shape2[0]/self.factor, nlines=shape2[1]/self.factor,
      #                interpolant=self.interpolant)
      xshift = shape2[0] / 2. - pm1[0]
      yshift = shape2[1] / 2. - pm1[1]
      output_data = pyfits.getdata(output)
      print "xshift, yshift:", xshift, yshift
      output2_data = interpolation.shift(output_data, 
                     [xshift, yshift], order=1, mode='wrap')
      # check the new psf is centered...
      pm2 = measurements.center_of_mass(output2_data)
      print "Shifted center of mass: %.2f, %.2f" % tuple(pm2)
      # Now zoom the PSF back to the original pixel scale
      output2_data = interpolation.zoom(output2_data, 1./self.factor, order=1,
                     mode='wrap')
      pyfits.append(output2, output2_data, pyfits.getheader(output))
      # Now center the PSF
      # pm2 = image_moments.moments(pyfits.getdata(output2))
      # npix = np.min(pyfits.getdata(output2).shape)
      # print "npix", npix
      # iraf.imlintran(input=output2, output='temp.fits',
      #                xrotation=0., yrotation=0.,
      #                xin=pm2.x1+1, yin=pm2.y1+1,
      #                xmag=1.0, ymag=1.0, 
      #                ncols=npix, nlines=npix, 
      #                interpolant='linear')
      # Now normalize the PSF
      # imgsum = pyfits.getdata(output2).ravel().sum()
      imgsum = output2_data.sum()
      imgsum = np.abs(imgsum)
      iraf.imarith(output2, '/', imgsum, 'temp.fits')
      os.remove(output2)
      os.system('mv temp.fits %s' % output2)
      os.system('rm temp*.fits')
      if not save:
         print "Removing %s..." % output
         os.remove(output)

   def run_all(self, mode='median'):
      self.cut_all()
      self.oversample_all()
      self.stack_all(mode=mode, save=self.save)
      print "Done."
