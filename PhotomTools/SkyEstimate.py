#!/usr/bin/env python

import numpy as np
import pyfits
from stats import gauss, robust
from ImageTools import images
import matplotlib.pyplot as plt


class SkyEstimate(object):
   """
   Estimate the sky noise directly from the image, instead of using RMS or 
   weight map. We fit a Gaussian to the sky pixel distribution, rejecting
   the positive tail from pixels contaminated by wings of nearby objects.
   """
   def __init__(self, drzimage, flagimage, filtername):
      self.drzimage = images.FITSImage(drzimage)
      self.flgimage = images.FITSImage(flagimage)
      self.filtername = filtername
      self.pixscale = self.drzimage.pixscale()  # in arcsec / pixel
      self.skyimage = np.where(self.flgimage.data==0, self.drzimage.data, 1.e10)

   def skyest_xy(self, x, y, width=64, minpix=100, nsigma=3.0, plot=True):
      xlo = int(round(x - width/2))
      xhi = int(round(x + width/2))
      ylo = int(round(y - width/2))
      yhi = int(round(y + width/2))
      # collect sky pixels
      skypix = self.skyimage[ylo:yhi,xlo:xhi]
      skypix = skypix[skypix<1.e10].ravel()
      if len(skypix) <= minpix:
         print "************************************************"
         print "Not enough sky pixels (%d) to determine noise..." % len(skypix)
         print "Return -1 for both background and RMS."
         print "************************************************"
         self.bkgd = -1.
         self.sigma = -1.
         return (-1., -1.)
      else:
         print "Number of sky pixels:", len(skypix)
         print "Start fitting Gaussian to sky pixel distribution..."
         bkgd, sigma = gauss.fitgauss(skypix, clip=True, clipsigma=nsigma)
         print "Estimated sky background level is %f" % bkgd
         print "Estimated sky noise level is %f" % sigma
         self.bkgd = bkgd
         self.sigma = sigma
         if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            bins = np.linspace(skypix.min(), skypix.max(), num=50)
            xarr = np.linspace(skypix.min(), skypix.max(), num=500)
            h = ax.hist(skypix, bins=bins, histtype='step', lw=2.0)
            g = gauss.gauss(xarr, bkgd, sigma)
            ax.plot(xarr, g * h[0].max() / g.max(), ls='--', lw=3.0)
            ax.set_title('Sky pixel around (%.2f, %.2f)' % (x, y), size=18)
            ax.set_xlabel('%s pixel value' % self.filtername.upper(), size=16)
            ax.set_ylabel('Number', size=16)
            skytxt = "bkgd = %.3e\nsigma = %.3e" % (bkgd, sigma)
            ax.text(0.95, 0.95, skytxt, ha='right', va='top', 
                    transform=ax.transAxes,
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
         return bkgd, sigma

   def skyest_radec(self, ra, dec, width=2.0, nsigma=3.0, plot=True):
      """
      Given the sky coordinate of an object, estimate the local sky noise 
      within a square box of width arcsec on a side, masking out pixels 
      belonging to the segmentations of detected objects.
      """
      # first calculate the (x, y) coordinate of the center
      # note that (x, y) should be 0-based (after python convention)
      x, y = self.drzimage.wcs.wcs_sky2pix([[ra, dec]], 0)[0]
      print "(x, y) = ", x, y
      width_pix = width / self.pixscale
      self.skyest_xy(x, y, width=width_pix, nsigma=nsigma, plot=plot)

   def calc_aperture_noise(self, apersize):
      """
      Calculate the total RMS within a circular aperture of diameter apersize
      (in arcsec).
      """
      rad_pix = (apersize / self.pixscale) / 2.0
      npix = 4. * np.pi * rad_pix**2
      # calculate total RMS by adding the variances from each pixel
      sigma_tot = np.sqrt(npix * self.sigma**2)
      return sigma_tot

   def skymap(self, backsize=64):
      """
      Calculate a sky map using Gaussian fitting. Estimate background value af
      a grid of points, and then interpolate over the whole image.
      """
      ncols, nlines = self.drzimage.data.shape
      nxgrid = int(round(ncols / backsize))
      nygrid = int(round(nlines / backsize))
      raise NotImplementedError

