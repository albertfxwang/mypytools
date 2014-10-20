#!/usr/bin/env python

import numpy as np
from pygoods import *
import os
import pyfits

"""
Given the TFIT output catalog *.cat_best and magnitude zero-points (along with
FLUXCONV keyword for IRAC images), calculate the magnitude & magnitude errors
for the TFIT output flux.
"""
fluxconv = {'ch1':0.1253, 'ch2':0.1469}
# the sky RMS in the primary field, taken from Marusa et al. 2013 (in prep)
skyrms = {'bullet':{'ch1':0.000931, 'ch2':0.00104},
          'macs1149':{'ch1':0.000905, 'ch2':0.00110},
          'macs0454':{'ch1':0.000992, 'ch2':0.00119},
          'macs0717':{'ch1':0.000957, 'ch2':0.00114},
          'macs0744':{'ch1':0.000904, 'ch2':0.00109},
          'macs1423':{'ch1':0.000757, 'ch2':0.00098},
          'macs2129':{'ch1':0.000840, 'ch2':0.00105},
          'macs2214':{'ch1':0.000890, 'ch2':0.00110}}

#def tfit_mag(catbest, cluster, band, cutoutdir, fluxunit='MJy', 
def tfit_mag(catbest, fluxunit='MJy',
             zeropoint=21.581):
   """
   Note that if the IRAC mosaic has a pixel scale different from the native
   scale of 1.2 arcsec, one needs to supply a value that is f * FLUXCONV, where
   f = (1.2 / pixscale) ** 2 and pixscale is the pixel scale of the IRAC 
   mosaic.
   Assume that the fluxes in catbest are measured from the image with unit of 
   DN/sec, so the total flux from catbest also has the unit of DN/sec. 
   The zeropoint supplied is the zeropoint when fluxes are in ** MJy/sr **, so 
   one should convert the fluxes in the TFIT catalog to MJy/sr first...
   If the best-fit flux is negative or zero, set its magnitude to 99.0 and use 
   its flux error as the 1-sigma upper limit (i.e., the value for magnitude 
   error is the 1-sigma magnitude upper limit).
   The default value of zeropoint is 21.58 mag, which is in the AB magnitude
   system (with F0 = 3631 Jy instead of 280.9 Jy in ch1 as used in IRAC 
   instrument handbook), for images with pixel scale of 0.6 arcsec/pixel.

   band: either 'ch1' or 'ch2'
   fluxunit: either 'MJy' or 'DNS' (which is DN/sec)
   bkgd_err: 1-sigma error of the median background; this number will be added
             to the TFIT-reported fitquerr in quadrature to calculate the total
             magnitude error.
   """
   c = sextractor(catbest)
   if fluxunit == 'DNS':
      flux_sb = c.fitqty * fluxconv[band]
   else:
      flux_sb = c.fitqty.copy()
   fluxerr = c.fitquerr.copy()
   ston = flux_sb / fluxerr
   if fluxunit == 'DNS':
      fluxerr = fluxerr * fluxconv[band]
   # convert the flux back to MJy/sr
   #mag = np.where(c.fitqty > 0, new_zpt - 2.5 * np.log10(c.fitqty), 99.0)
   # mag = np.where(flux_sb > 0, zeropoint - 2.5 * np.log10(flux_sb), 99.0)
   mag = np.where(ston>1.0, zeropoint - 2.5 * np.log10(flux_sb), 99.0)
   # ston = c.fitqty / fluxerr
   #mag_err = np.where(ston > 0, 1.0857 / ston, 
   #                   new_zpt - 2.5 * np.log10(c.fitquerr))
   mag_err = np.where(ston > 1, 2.5*np.log10(1+1./ston),
                      zeropoint - 2.5 * np.log10(fluxerr))
   # now write to an output
   f = open(catbest+'_mag', 'wb')
   f.write(c._header)
   ncols = len(c._colnames)
   f.write('# %d mag\n' % (ncols+1))
   f.write('# %d magerr\n' % (ncols+2))
   #f.write('# %d fluxerr_tot\n' % (ncols+3))
   #f.write('# %d npix_cutout\n' % (ncols+4))
   for j in range(len(c)):
      for i in range(ncols):
         fmtstr = c._fmt[c._colnames[i]]     
         f.write(fmtstr % c.__getattribute__(c._colnames[i])[j])
         f.write(' ')
      # Now write magnitudes and magnitude errors
      f.write('%.4f %.4f ' % (mag[j], mag_err[j]))
      #f.write('%.6f ' % fluxerr[j])
      #f.write('%d ' % npix_cut[j])
      f.write('\n')

   f.close()