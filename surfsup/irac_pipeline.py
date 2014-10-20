#!/usr/bin/env python

import numpy as np
import pyfits
from pyraf import iraf
iraf.stsdas()
import subtract_scattered_background as ssb
import rmscheck
import os, sys
import pylab as plt
from stats import gauss

"""
A pipeline that prepares the IRAC mosaics for use by TFIT. 
Assume the starting images to be the following:
- a drizzled IRAC mosaic image with pixel scale 0.6 arcsec/pixel, with the 
  same astrometry as the high-resolution HST mosaics, and has the unit of 
  MJy/sr
- an IRAC uncertainty image (*_unc.fits), in units of MJy/sr
- usually the cutout of the IRAC mosaic should cover a larger area than the HST
  mosaics, so that TFIT doesn't break down when there are high-resolution 
  templates near (or even beyond) the image border of the IRAC mosaic
"""
exptimes={'ch1':96.8, 'ch2':96.8, 'ch3':96.8, 'ch4':46.8}   # exposure time per AOR (used to convert cov map to exp map)
fluxconv = {'ch1':0.1253, 'ch2':0.1469}
subreg_default = {'bullet':(141,260,365,450)}
# the FLUXCONV factor in units of MJy/sr per DN/sec, for the warm mission only
# taken from the webpage 
# http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/warmfeatures/

def irac_pipeline(drz_img, unc_img, channel, rmsmode='unc', pixscale=0.6, gain=3.7,
                  zeromask=None, bgmask=None, source_threshold=3,
                  mask_growsig=5, subreg='none', aperture_width=7,
                  nnear=23, constant_bkgd=False, cluster='bullet',
                  inmask_checkrms='ch1_pf0.01_cut_mask.fits'):
   """
   drz_img: the file name of the IRAC science mosaic (in MJy/sr)
   unc_img: the file name of the IRAC uncertainty image (in MJy/sr)
   channel: the IRAC channel name (either 'ch1' or 'ch2')
   pixscale: the pixel scale of the input IRAC mosaic in arcsec/pixel
   gain: the CCD gain of the IRAC mosaic
   zeromask: the image name that provides which pixels to be zeroed out at the 
             end of this pipeline, both in the drz image and the RMS image
   bgmask: user-provided mask for background subtraction (optional)
   subreg: provide a list (or array) of four numbers (xmin, xmax, ymin, ymax)
           that specifies a region in which to measure RMS. If None, default
           to the entire image?
   rmsmode: either 'cov' (use the coverage map to calculate the RMS map) or 
            'unc' (use the uncertainty map from the MOPEX pipeline)
   """
   clean_pipeline(drz_img)
   print "Step 1: convert image units from MJy/sr to DN/sec..."   
   drz_new1 = os.path.split(drz_img)[-1][:-5] + '_MJy.fits'  # without converting to DN/sec
   iraf.imcopy(drz_img, drz_new1)

   if zeromask != None:
      # zero out the masked pixel before entering background subtraction stage
      iraf.imcalc("%s,%s" % (drz_new1, zeromask), "temp1.fits", 
                  "if im2 > 0 then 0. else im1")
      os.remove(drz_new1)
      os.system('mv temp1.fits %s' % drz_new1)

   # Run background-subtraction step
   print "Step 2: subtract background from the science mosaic..."
   print "See ch1_sub_bkgd.param for explanations of each parameter."
   # Use the default parameters here... can maybe fiddle with these numbers 
   # later?
   drz_root2 = os.getcwd() + '/' + os.path.split(drz_new1)[-1][:-5]
   drz_new3 = drz_root2 + '_bgcorr.fits'
   if constant_bkgd == True:
      print "Subtracting a constant background within the designated region..."
      if bgmask == None:
         print "Error: must supply a background mask if one wants to subtract a constant background..."
         sys.exit()
      else:
         subtract_constant_bkgd(drz_new1, bgmask, drz_new3, 
                                growsig=mask_growsig)
   else:
      bkgd_parfile = os.path.split(drz_new1)[-1][:-5] + '_bkgd.param'
      print "bkgd_parfile", bkgd_parfile
      f = open(os.getcwd() + '/' + bkgd_parfile, 'wb')
      #drz_root2 = drz_new1[:-5]
      
      f.write('IMAGE \t %s\n' % drz_new1)
      f.write('OUT_MASK \t %s\n' % (drz_root2 + '_bgmask.fits'))
      f.write('OUT_BKGD \t %s\n' % (drz_root2 + '_bkgd.fits'))
      f.write('OUT_BKGD_SUBTRACTED \t %s\n' % (drz_root2 + '_bgcorr.fits'))
      f.write('OUT_NEWCOLLAGE \t %s\n' % (drz_root2 + '_bkgdcoll.fits'))
      f.write('NNEAR \t %d\n' % nnear)
      f.write('SOURCE_THRESHOLD \t %f\n' % float(source_threshold))
      f.write('MASK_GROWSIG \t %f\n' % float(mask_growsig))
      f.write('DQEXT \t 0\n')
      f.write('CRUDESUB \t 1\n')
      if subreg == None:
         f.write('SUBREG \t %d,%d,%d,%d\n' % tuple(subreg_default[cluster]))
      elif subreg == "none":
         f.write('SUBREG \t none\n')
      else:
         f.write('SUBREG \t %d,%d,%d,%d' % tuple(subreg))
      f.write('\n')
      f.write('DILATE_FACTOR \t 2\n')
      f.write('DILATE_THRESH \t 300\n')
      f.write('DILATE_TOT_THRESH \t 500\n')
      f.write('DILATE_SATURATED \t 1\n')
      f.write('GAIN \t %f\n' % gain)
      f.write('OUT_LABEL \t %s\n' % (drz_root2 + '_bgmask_labels.fits'))
      f.write('\n')
      f.write('MAKE_RMS \t 0\n')
      f.write('OUT_RMS \t %s\n' % (drz_root2 + '_rms_estimated.fits'))
      f.write('RMS_SMOOTH \t 100\n')   
      f.write('\n')
      f.write('OUT_CRUDE \t %s\n' % (drz_root2 + '_simple_bgcorr.fits'))
      f.write('OUT_FILT  \t %s\n' % (drz_root2 + '_bkgd.filt.fits'))
      if bgmask != None:
         f.write('IN_MASK \t %s\n' % bgmask)
      f.write('APERTURE_WIDTH \t %d\n' % aperture_width)
      f.write('\n')
      f.close()
      ssb.run(bkgd_parfile)
      # Make an image of background pixels after subtraction
      #if os.path.exists(drz_root2+'_bgpix.fits'):
      #   os.remove(drz_root2+'_bgpix.fits')
      #iraf.imcalc('%s,%s' % (drz_root2+'_bgcorr.fits','ch1_pf0.01_cut_seg.fits'),
      #            drz_root2+'_bgpix.fits','if im2 > 0 then -99.0 else im1')
   
   if zeromask != None:
      # again, zero out the masked pixels in the background-subtracted image
      print "zero out the masked pixels in the background-subtracted image..."
      iraf.imcalc("%s,%s" % (drz_new3, zeromask), "temp2.fits", 
                  "if im2 > 0 then 0. else im1")
      os.remove(drz_new3)
      os.system('mv temp2.fits %s' % drz_new3)


   # print "Step 3: re-scale the uncertainty image"
   if rmsmode == 'cov':
      print "Step 3: make RMS map from coverage map, then re-scale the uncertainty image"
      cov_img = unc_img
      rms_new1 = drz_img[:-9] + '_rms.fits'
      if os.path.exists(rms_new1):
         os.remove(rms_new1)
      iraf.imcalc(cov_img, rms_new1, 
                  "if im1 == 0. then 1000. else 1./sqrt(im1*%f)" % exptimes[channel])
      # use coverage map (in units of seconds) as the weight map to calculate
      # the RMS map
      unc_new1 = rms_new1  # to comply with previous version
   else:
      print "Step 3: use the uncertainty map as the RMS map, and use FLUXCONV to convert to DN/sec "
      print "(the important point is to convert the RMS map to the same unit as the science map)"
      #unc_new1 = unc_img[:-5] + '_DNS.fits'
      unc_new1 = unc_img[:-5] + '_MJy.fits'
      #iraf.imcalc(unc_img, unc_new1, 'im1/%.4f' % (factor * fluxconv[channel]))
      #iraf.imcalc(unc_img, unc_new1, 'im1/%.4f' % (fluxconv[channel]))
      iraf.imcopy(unc_img, unc_new1)  # without converting to DN/sec
   try:
      inmask = pyfits.getdata(inmask_checkrms)
   except:
      inmask = None
   rmscheck.checkrms(pyfits.getdata(drz_new3), pyfits.getdata(unc_new1),
                     inmask=inmask, growsig=0.2)
   try:
      unc_scale = float(raw_input('Now enter the scaling factor for the uncertainty image: '))
   except:
      print 'Need to enter a number here.'
      sys.exit()
   unc_new2 = unc_new1[:-5] + '_scaled.fits'
   if os.path.exists(unc_new2):
      os.remove(unc_new2)
   print "Multiply the uncertainty image by a factor of %f" % unc_scale
   iraf.imcalc(unc_new1, unc_new2, 'im1*%f' % unc_scale)

   print "Step 4: make the source-weighted RMS map..."
   # Here both the science image and the (re-scaled) RMS map are in units of 
   # DN/sec
   #mathstr = 'sqrt(im1**2 + im2/%f)' % gain
   #mathstr = 'sqrt(im1**2 + im2)'  # a different weighting scheme
   #mathstr = 'sqrt((im1*%f)**2 + im2*%f) / %f' % (gain, gain, gain)   
   #mathstr = 'sqrt((im1*%f)**2 + im2) / %f' % (gain, gain)
   mathstr = 'im1 * 1.0'
   # convert everything into photon units
   unc_new4 = unc_new2[:-5] + '_src.fits'
   if os.path.exists(unc_new4):
      os.remove(unc_new4)
   iraf.imcalc('%s,%s' % (unc_new2,drz_new3), unc_new4, mathstr)
   #iraf.imreplace(unc_new4, 1000., lo=0., up=0.)
   iraf.imcalc('%s,%s' % (unc_new4,unc_new2), 'temp_unc.fits', 
              'if im1 <= 0. then im2 else im1')
   os.system('mv temp_unc.fits %s' % unc_new4)   
   # Extra steps


   print "All steps finished. A summary:"
   print "Background-subtracted science mosaic: %s" % drz_new3
   print "Source-weighted, re-scaled RMS map: %s" % unc_new4
   print "All images converted into DN/sec."

def clean_pipeline(drzimg, suffix='DNS'):
   """
   Clean the products of IRAC pipeline, given the input mosaic.
   """
   root1 = os.path.splitext(drzimg)[0]
   root2 = root1[:-4]  # assume that root ends with "_drz"
   os.system('rm %s_drz_%s*.fits' % (root2, suffix))
   os.system('rm %s_unc_%s*.fits' % (root2, suffix))
   os.system('rm %s_%s_bkgd.param' % (root1, suffix))

def subtract_constant_bkgd(input_image, mask_image, output_image, growsig=0,
                           subreg=None, savemask=True, boxsize=40, 
                           boxcenter=None, fitgauss=True, show_plot=False):
   """
   Subtract a constant (median) background not masked by the mask_image, where 
   mask_image > 0 means the pixels are masked. The mask image should contain
   at least the segmentations of sources, and could also define the subregion
   of an image.

   boxsize: the size of the square box centered around the object within which
            one determines the median background
   """
   #if subreg == None:
   #   img = pyfits.getdata(input_image)
   #   subreg = [1, img.shape[0], 1, img.shape[1]]
   # output_image = os.getcwd() + '/' + os.path.split(output_image)[-1]
   img1 = pyfits.getdata(input_image)
   mask = pyfits.getdata(mask_image)
   if growsig > 0:
      mask = ssb.growmask(mask, growsig=growsig)
      mask = np.where(mask, 1, 0)
      maskout = os.path.splitext(mask_image)[0] + '_out.fits'
      if os.path.exists(maskout):
         os.remove(maskout)
      hdu = pyfits.PrimaryHDU(mask)
      hdu.writeto(maskout)
   if subreg != None:
      xmin, xmax, ymin, ymax = subreg
   elif boxsize != None:
      xc, yc = boxcenter
      xmin = xc - boxsize / 2
      xmax = xc + boxsize / 2
      ymin = yc - boxsize / 2
      ymax = yc + boxsize / 2
   else:
      xmin = 0
      xmax = img1.shape[0]
      ymin = 0
      ymax = img1.shape[1]
   print "Determining median background within the region:"
   print "(xmin, xmax, ymin, ymax) = (%d, %d, %d, %d)" % (xmin, xmax, ymin, ymax)
   mask_subreg = np.ones(img1.shape, 'int')
   mask_subreg[ymin-1:ymax,xmin-1:xmax] = 0
   mask = np.logical_or(mask, mask_subreg)
   mask = np.where(mask==True, 1, 0)
   img1_masked = np.ma.masked_array(img1, mask=mask)
   bgpix = np.compress(mask.ravel()==0, img1.ravel())
   y = np.histogram(bgpix, bins=30)
   # figure out if the best-fit peak is too far away from the peak of the 
   # histogram, which suggests bimodality
   ipeak = np.argsort(y[0])[-1]
   # calculate the bin center of the peak of histogram
   xmode = y[1][ipeak] + (y[1][1] - y[1][0]) / 2.
   if fitgauss:
      ## Determin local background by fitting a Gaussian...
      bkgd_median, bkgd_sig = gauss.fitgauss(bgpix, clip=True, clipsigma=3.0)
      # if np.abs(bkgd_median - xmode) >= 0.5 * bkgd_sig:
      # do something...?
   else:
      ## ... or by just calculating the median & standard deviation
      bkgd_median = np.ma.median(img1_masked)
      bkgd_sig = np.ma.std(img1_masked)
   num_bkgd = np.sum(mask == 0)
   print "Median sky background: %.6f" % bkgd_median
   print "1-sigma of sky background: %.6f" % bkgd_sig
   print "Number of sky pixels for determining background: %d" % num_bkgd
   if os.path.exists(output_image):
      os.remove(output_image)
   iraf.imcalc(input_image, output_image, 'im1 - %f' % bkgd_median)
   h = pyfits.open(output_image, mode='update')
   h[0].header['MEDBKGD'] = bkgd_median
   h[0].header['SIGBKGD'] = bkgd_sig
   h[0].header['NUMBKGD'] = num_bkgd
   h.flush()
   h.close()
   # Show background pixel histogram
   print "Number of background pixels:", len(bgpix)
   fig = plt.figure()
   ax = fig.add_subplot(111)
   y = ax.hist(bgpix, bins=y[1], histtype='step', lw=2.0)
   # also plot the best-fit gaussian
   x = np.linspace(bgpix.min(), bgpix.max(), num=100)
   g = gauss.gauss(x, bkgd_median, bkgd_sig)
   g = g * np.max(y[0]) / g.max()
   ax.plot(x, g, color='red', lw=2.0)
   bkgd_str = "median=%.2e\nsigma=%.2e" % (bkgd_median, bkgd_sig)
   ax.text(0.95, 0.95, bkgd_str, transform=ax.transAxes, ha='right', va='top',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
   ymax = ax.get_ylim()[1]
   ax.plot([bkgd_median, bkgd_median], [0., ymax], ls='--', lw=2.0, color='black')
   plt.title('Background pixel values', size=24)
   if not show_plot:
      plt.close('all')
   if savemask:
      maskname = os.path.splitext(output_image)[0] + '_mask.fits'
      if os.path.exists(maskname):
         os.remove(maskname)
      hdu2 = pyfits.PrimaryHDU(mask)
      hdu2.writeto(maskname)
   return fig

