#!/usr/bin/env python

import numpy as np
from tfit_tools import align_images
import yaml
from pyraf import iraf
iraf.images()
import pyfits
import os, sys, glob
import subprocess
import interp_iracpsf
import make_tfitcat
from ImageTools import images

"""
A pipeline that prepares the HST and IRAC images for running TFIT. An example
parameter file is at $mypytools/surfsup/example_tfit_prep.yml

Before running this pipeline, make cutouts of the HST images to trim the borders.
Do the same for IRAC images, and make sure that IRAC coverage is larger than
HST coverage.

Step 1: Re-sample the HST images so that the pixel scales between HST and IRAC
images are integer multiples. Also take care of CRVALs and half-integer CRPIXs.

Step 2: Run SExtractor on HST images to prepare a high-resolution catalog. 
Runs background subtraction on HST images (optional).

Step 3: Run SExtractor on IRAC images to make a segmentation map, which will 
serve as a mask when calculating local background for IRAC sources.
"""
default_bkgd_params = '/Users/khuang/Dropbox/Research/surfsup_dropbox/example_bkgd.param.yml'

class TFITPrep(object):
   def __init__(self, paramfile):
      c = yaml.load(open(paramfile))
      self.cluster_name = c['cluster_name']
      self.hr_dir = c['hr_dir']
      self.lr_dir = c['lr_dir']
      self.stages = c['stages']
      # self.swarp_params = []
      self.hr_bands = c['hr_bands']
      self.lr_bands = c['lr_bands']
      if 'swarp' in self.stages:
         self.swarp_file = c['swarp_file']
      self.hr_scale = c['hr_scale']
      self.newscale = '%2d' % int(c['hr_scale']*1000.)
      self.root = '%s_4tfit' % (self.cluster_name)
      # self.drzroot = os.path.splitext(c['hst_photom_drz'])[0]
      # Setting up file names for hires images
      self.hr_input_drz = dict(zip(self.hr_bands, c['hr_input_drz']))
      self.hr_input_wht = dict(zip(self.hr_bands, c['hr_input_wht']))
      self.hr_input_flg = c['hr_input_flg']
      self.hr_output_drz = dict(zip(self.hr_bands, c['hr_output_drz']))
      self.hr_output_wht = dict(zip(self.hr_bands, c['hr_output_wht']))
      self.hr_output_flg = c['hr_output_flg']
      self.hr_output_seg = c['hr_output_seg']
      if 'blocksum_factor' in c.keys():
         self.blocksum_factor = c['blocksum_factor']
      else:
         self.blocksum_factor = 1
      # self.hr_resample_drz = dict(zip(self.hr_bands, 
      #    map(lambda b: '%s_%s_drz.fits' % (b, self.root), self.hr_bands)))
      # self.hr_resample_wht = dict(zip(self.hr_bands, 
      #    map(lambda b: '%s_%s_wht.fits' % (b, self.root), self.hr_bands)))
      # Flag and seg images only for the photometry band
      # self.hr_resample_flg = '%s_%s_flg.fits' % (self.hr_bands[-1].upper(), self.root)
      # self.hr_resample_seg = '%s_%s_seg.fits' % (self.hr_bands[-1].upper(), self.root)
      # else:
      #    self.hr_resample_drz = self.hr_input_drz
      #    self.hr_resample_wht = self.hr_input_wht
      #    self.hr_resample_flg = c['hr_input_flg']
      #    self.hr_resample_seg = c['hr_input_seg']
      self.hr_sexcat = c['hr_sexcat']
      self.hr_tfit_cat = c['hr_tfit_cat']
      # Set up lo-res image names
      self.lr_drz = dict(zip(self.lr_bands, c['lr_drz']))
      self.lr_unc = dict(zip(self.lr_bands, c['lr_unc']))
      self.lr_wht = dict(zip(self.lr_bands, c['lr_wht']))
      self.lr_flg = dict(zip(self.lr_bands, c['lr_flg']))
      self.use_bgcorr = c['use_bgcorr']
      if len(self.hr_bands) > 1:
         self.bgcorr_dir = '%s_bgcorr' % self.hr_bands[1]
         self.bgcorr_image = '%s/%s_%s_bgcorr.fits' % (self.bgcorr_dir, self.hr_bands[1], self.root)
      else:
         self.bgcorr_dir = '%s_bgcorr' % self.hr_bands[0]
         self.bgcorr_image = '%s/%s_%s_bgcorr.fits' % (self.bgcorr_dir, self.hr_bands[0], self.root)
      self.magzpt = c['magzpt']
      self.hr_sexfile = c['hr_sexfile']
      self.lr_sexfile = c['lr_sexfile']
      os.chdir(self.hr_dir)


   def swarp(self):
      ### Runs Swarp on high-res images; can be run on more than 2 bands.
      # for sp in self.swarp_params:
      #    align_images.swarp_images(sp)
      for b in self.hr_bands:
         sp = {}
         sp['swarp_file'] = self.swarp_file
         sp['hires_dir'] = self.hr_dir
         # sp['lores_dir'] = os.path.join(self.lr_dir, self.lr_bands[0])
         sp['lores_dir'] = self.lr_dir
         sp['hires_input_drz'] = self.hr_input_drz[b]
         sp['hires_input_wht'] = self.hr_input_wht[b]
         # sp['hires_output_drz'] = self.hr_resample_drz[b]
         # sp['hires_output_wht'] = self.hr_resample_wht[b]
         sp['hires_output_drz'] = self.hr_output_drz[b]
         sp['hires_output_wht'] = self.hr_output_wht[b]
         sp['hr_scale'] = self.hr_scale
         sp['lores_drz'] = self.lr_drz[self.lr_bands[0]]
         sp['lores_unc'] = self.lr_unc[self.lr_bands[0]]
         sw_name = '%s_%s.swarp.yml' % (b, self.cluster_name.lower())
         yaml.dump(sp, open(sw_name, 'wb'), default_flow_style=False)
         align_images.swarp_images(sw_name)
         # Make flag image for photometry band
         if b == self.hr_bands[1]:
            ### Calculate flag image for the detection band (=self.hr_bands[1])
            if os.path.exists(self.hr_output_flg):
               os.remove(self.hr_output_flg)
            # wht_img = pyfits.getdata(self.hr_resample_wht[b]).astype(np.float)
            wht_img = pyfits.getdata(self.hr_output_wht[b]).astype(np.float)
            print "wht_img.dtype", wht_img.dtype
            # wht_hdr = pyfits.getheader(self.hr_resample_wht[b])
            wht_hdr = pyfits.getheader(self.hr_output_wht[b])
            flg_img = np.where(wht_img > 0, 0, 1).astype(np.int32)
            print "flg_img.dtype", flg_img.dtype
            # pyfits.append(self.hr_resample_flg, flg_img, wht_hdr)
            pyfits.append(self.hr_output_flg, flg_img, wht_hdr)
            # flg_hdu = pyfits.PrimaryHDU(flg_img)
            # for k in wht_hdr.keys():
            #    if k not in flg_hdu.header.keys():
            #       try:
            #          flg_hdu.header.set(k, wht_hdr[k])
            #       except:
            #          print "Unable to add keyword %s to the header; skip..." % k
            # flg_hdulist = pyfits.HDUList([flg_hdu])
            # flg_hdulist.writeto(self.hr_resample_flg)
      # Now shift the WCS for the other low-res bands as well
      lr_hdr = pyfits.getheader(os.path.join(self.lr_dir, self.lr_drz[self.lr_bands[0]]))
      crval1 = lr_hdr['crval1']
      crval2 = lr_hdr['crval2']
      crpix1 = lr_hdr['crpix1']
      crpix2 = lr_hdr['crpix2']
      crvals_new = [crval1, crval2]
      crpixs_new = [crpix1, crpix2]
      for lb in self.lr_bands[1:]:
         align_images.update_crvals_crpixs(os.path.join(self.lr_dir, self.lr_drz[lb]), crvals_new, crpixs_new)
         align_images.update_crvals_crpixs(os.path.join(self.lr_dir, self.lr_unc[lb]), crvals_new, crpixs_new)

   def blocksum(self):
      ### If the input images already have an integer multiple pixel scale from
      ### the IRAC mosaics, but one wants to create larger pixels to aid 
      ### faint-source detections...
      #raise NotImplementedError
      # For science images, just block sum
      print "Block-summing the science images..."
      for b in self.hr_bands:
         print b.upper()
         # iraf.blkavg(self.hr_input_drz[b], self.hr_output_drz[b], 
         #             self.blocksum_factor, self.blocksum_factor, option='sum')
         drz = images.FITSImage(self.hr_input_drz[b])
         drz.blocksum(self.blocksum_factor)
         drz.writeto(self.hr_output_drz[b], overwrite=True)
      # Now sum the variance images, and then convert back to weight images
      print "Block-summing the weight images (actually the variance image)..."
      for b in self.hr_bands:
         print b.upper()
         hr_input_var = self.hr_input_wht[b].replace('wht','var')
         if os.path.exists(hr_input_var):
            os.remove(hr_input_var)
         # Make a variance image; first make a copy of the weight image
         os.system('cp %s %s ' % (self.hr_input_wht[b], hr_input_var))
         # then calculate var = 1./wht
         h = pyfits.open(hr_input_var, mode='update')
         h[0].data = np.where(h[0].data > 0, 1. / h[0].data, 1.e10)
         h.flush()
         h.close()
         # Now block sum the variance image; use the same image name so the 
         # block-summed variance image replace the old variance image
         print hr_input_var
         # iraf.blkavg(hr_input_var, hr_input_var, 
         #             self.blocksum_factor, self.blocksum_factor, option='sum')
         var = images.FITSImage(hr_input_var)
         var.blocksum(self.blocksum_factor)
         var.writeto(name=hr_input_var, overwrite=True)
         # Now convert back to weight image
         os.system('cp %s %s ' % (hr_input_var, self.hr_output_wht[b]))
         h = pyfits.open(self.hr_output_wht[b], mode='update')
         h[0].data = np.where(h[0].data < 1.e10, 1. / h[0].data, 0.)
         h.flush()
         h.close()
      # Make flag images
      print "Make flag images from weight images..."
      for b in self.hr_bands:
         print b.upper()
         wht = images.WeightImage(self.hr_output_wht[b])
         wht.make_flg(self.hr_output_wht[b].replace('wht', 'flg'))


   def shift_crpix(self):
      hb0 = self.hr_bands[0]
      if ('swarp' not in self.stages) and ('blocksum' not in self.stages):
         print "Make copies of the input images..."
         if self.hr_input_drz[hb0] != self.hr_output_drz[hb0]:
            os.system('cp %s %s ' % (self.hr_input_drz[hb0], self.hr_output_drz[hb0]))
         if self.hr_input_wht[hb0] != self.hr_output_wht[hb0]:
            os.system('cp %s %s ' % (self.hr_input_wht[hb0], self.hr_output_wht[hb0]))
         if self.hr_input_flg != self.hr_output_flg:
            os.system('cp %s %s ' % (self.hr_input_flg, self.hr_output_flg))
         # os.system('cp %s %s ' % (self.hr_input_drz[hb0], self.hr_resample_drz[hb0]))
         # os.system('cp %s %s ' % (self.hr_input_wht[hb0], self.hr_resample_wht[hb0]))
         # os.system('cp %s %s ' % (self.hr_input_flg, self.hr_resample_flg))
      print "Shifting CRPIXs for %s..." % hb0
      for lb in self.lr_bands:
         align_images.shift_crpix(self.hr_output_drz[hb0], 
                                  self.hr_output_wht[hb0],
                                  self.lr_drz[lb], self.lr_wht[lb], 
                                  hr_dir=self.hr_dir, lr_dir=self.lr_dir)
         # align_images.shift_crpix(self.hr_resample_drz[hb0], 
         #                          self.hr_resample_wht[hb0],
         #                          self.lr_drz[lb], self.lr_wht[lb], 
         #                          hr_dir=self.hr_dir, lr_dir=self.lr_dir)
      if len(self.hr_bands) > 1:
         for hb in self.hr_bands[1:]:
            # shift the CRVALs and CRPIXs to match the detection band
            # hb2 = self.hr_bands[1]
            print "Shifting CRPIXs for %s..." % hb
            h1 = pyfits.getheader(self.hr_output_drz[hb0])
            h2 = pyfits.open(self.hr_output_drz[hb], mode='update')
            h2[0].header['crval1'] = h1['crval1']
            h2[0].header['crval2'] = h1['crval2']
            h2[0].header['crpix1'] = h1['crpix1']
            h2[0].header['crpix2'] = h1['crpix2']
            h2.flush()
            h2.close()

   def hst_bkgd_sub(self, **newpars):
      """
      If values other than the default ones are desired, provide them as the 
      additional keyword arguments.
      """
      ### BROKEN FOR NOW... (5/19/2014)
      import subtract_scattered_background as ssb
      # Background subtraction for photometry band only
      bpars = yaml.load(open(default_bkgd_params))
      if not os.path.isdir(self.bgcorr_dir):
         os.mkdir(self.bgcorr_dir)
      if len(glob.glob('%s/*.fits' % self.bgcorr_dir)):
         os.system('rm %s/*.fits' % self.bgcorr_dir)
      # bpars['IMAGE'] = self.hr_resample_drz[self.hr_bands[-1]]
      bpars['IMAGE'] = self.hr_output_drz[self.hr_bands[-1]]
      bpars['OUT_MASK'] = '%s/%s_%s_bgmask.fits' % (self.bgcorr_dir, self.hr_bands[-1], self.root)
      bpars['OUT_BKGD'] = '%s/%s_%s_bkgd.fits' % (self.bgcorr_dir, self.hr_bands[-1], self.root)
      bpars['OUT_BKGD_SUBTRACTED'] = self.bgcorr_image
      bpars['OUT_NEWCOLLAGE'] = '%s/%s_%s_bkgdcoll.fits' % (self.bgcorr_dir, self.hr_bands[-1], self.root)
      bpars['OUT_LABEL'] = '%s/%s_%s_bgmask_labels.fits' % (self.bgcorr_dir, self.hr_bands[-1], self.root)
      if len(newpars.keys()):
         for k in newpars.keys():
            bpars[k.upper()] = str(newpars[k])
      self.bkgd_param = '%s_%s_bkgd.param' % (self.hr_bands[-1], self.root)
      f = open(self.bkgd_param, 'wb')
      strfmt = '%-24s %-29s\n'
      f.write(strfmt % ('IMAGE', bpars['IMAGE']))
      f.write(strfmt % ('OUT_MASK', bpars['OUT_MASK']))
      f.write(strfmt % ('OUT_BKGD', bpars['OUT_BKGD']))
      f.write(strfmt % ('OUT_BKGD_SUBTRACTED', bpars['OUT_BKGD_SUBTRACTED']))
      f.write(strfmt % ('OUT_NEWCOLLAGE', bpars['OUT_NEWCOLLAGE']))
      f.write(strfmt % ('NNEAR', bpars['NNEAR']))
      f.write(strfmt % ('SOURCE_THRESHOLD', bpars['SOURCE_THRESHOLD']))
      f.write(strfmt % ('MASK_GROWSIG', bpars['MASK_GROWSIG']))
      f.write(strfmt % ('DQEXT', bpars['DQEXT']))
      f.write(strfmt % ('CRUDESUB', bpars['CRUDESUB']))
      if type(bpars['SUBREG']) == type([]):
         f.write(strfmt % ('SUBREG', ','.join(bpars['SUBREG'])))
      else:
         f.write(strfmt % ('SUBREG', bpars['SUBREG']))
      f.write(strfmt % ('DILATE_FACTOR', bpars['DILATE_FACTOR']))
      f.write(strfmt % ('DILATE_THRESH', bpars['DILATE_THRESH']))
      f.write(strfmt % ('DILATE_TOT_THRESH', bpars['DILATE_TOT_THRESH']))
      f.write(strfmt % ('DILATE_SATURATED', bpars['DILATE_SATURATED']))
      f.write(strfmt % ('GAIN', bpars['GAIN']))
      f.write(strfmt % ('OUT_LABEL', bpars['OUT_LABEL']))
      f.write(strfmt % ('MAKE_RMS', bpars['MAKE_RMS']))
      f.write(strfmt % ('OUT_RMS', bpars['OUT_RMS']))
      f.write(strfmt % ('RMS_SMOOTH', bpars['RMS_SMOOTH']))
      f.close()
      ssb.run(self.bkgd_param)
      # Now zero-out the background-subtraction patterns outside the footprint
      h_bgcorr = pyfits.open(bpars['OUT_BKGD_SUBTRACTED'], mode='update')
      img_bgcorr = h_bgcorr[0].data
      # h_wht = pyfits.getdata(self.hr_resample_wht[self.hr_bands[-1]])
      h_wht = pyfits.getdata(self.hr_output_wht[self.hr_bands[-1]])
      img_bgcorr_zeroed = np.where(h_wht>0., img_bgcorr, 0.)
      h_bgcorr[0].data = img_bgcorr_zeroed
      h_bgcorr.flush()
      h_bgcorr.close()

   def hst_sextractor(self, **extra_kwargs):
      """
      Only run SExtractor on the photometry band in HST.
      """
      os.chdir(self.hr_dir)
      if self.use_bgcorr:
         input_image = self.bgcorr_image
      else:
         # input_image = self.hr_resample_drz[self.hr_bands[-1]]
         input_image = self.hr_output_drz[self.hr_bands[1]]
      sex_cmd = ['cex', '%s,%s' % (self.hr_output_drz[self.hr_bands[0]], input_image)]
      # sex_cmd = ['cex', '%s,%s' % (self.hr_resample_drz[self.hr_bands[0]], input_image)]
      sex_cmd += ['-c', self.hr_sexfile]
      sex_cmd += ['-WEIGHT_TYPE', 'MAP_WEIGHT,MAP_WEIGHT']
      sex_cmd += ['-WEIGHT_IMAGE', '%s,%s' % (self.hr_output_wht[self.hr_bands[0]], self.hr_output_wht[self.hr_bands[1]])]
      # sex_cmd += ['-WEIGHT_IMAGE', '%s,%s' % (self.hr_resample_wht[self.hr_bands[0]], self.hr_resample_wht[self.hr_bands[1]])]
      sex_cmd += ['-WEIGHT_THRESH', '0.']
      sex_cmd += ['-CATALOG_NAME', self.hr_sexcat]
      sex_cmd += ['-FLAG_IMAGE', self.hr_output_flg]
      # sex_cmd += ['-FLAG_IMAGE', self.hr_resample_flg]
      sex_cmd += ['-CHECKIMAGE_NAME', self.hr_output_seg]
      # sex_cmd += ['-CHECKIMAGE_NAME', self.hr_resample_seg]
      # sex_cmd += ['-PIXEL_SCALE', str(self.pixscale)]
      sex_cmd += ['-MAG_ZEROPOINT', str(self.magzpt)]
      if len(extra_kwargs):
         for k in extra_kwargs.keys():
            sex_cmd += ['-%s' % k.upper(), str(extra_kwargs[k])]
      print sex_cmd
      subprocess.call(sex_cmd)
      # Now format the SExtractor catalog into the TFIT-accepted form
      self.make_tfitcat()

   def make_tfitcat(self, hr_sexcat=None):
      if hr_sexcat == None:
         hr_sexcat = self.hr_sexcat
      make_tfitcat.make_tfitcat(hr_sexcat, outname=self.hr_tfit_cat)

   def irac_sextractor(self, **extra_kwargs):
      for b in self.lr_bands:
         wkdir = os.path.join(self.lr_dir, b)
         os.chdir(wkdir)
         print "curdir:", os.getcwd()
         drz_b = os.path.split(self.lr_drz[b])[-1]
         wht_b = os.path.split(self.lr_wht[b])[-1]
         flg_b = os.path.split(self.lr_flg[b])[-1]
         unc_b = os.path.split(self.lr_unc[b])[-1]
         # First, make weight images if they do not already exist
         if not os.path.exists(wht_b):
            unc_img = pyfits.getdata(unc_b)
            unc_hdr = pyfits.getheader(unc_b)
            wht_img = np.where(unc_img > 0., 1./unc_img**2, 0.).astype('float64')
            pyfits.append(wht_b, wht_img, unc_hdr)
         # Make flag image if not already there
         if not os.path.exists(flg_b):
            wht_img = pyfits.getdata(wht_b)
            wht_hdr = pyfits.getheader(wht_b)
            flg_img = np.where(wht_img > 0, 0, 1).astype('int32')
            pyfits.append(flg_b, flg_img, wht_hdr)
         broot = '%s_%s' % (b, self.cluster_name)
         # segimg = broot + '_seg.fits'
         segimg = wht_b.replace('wht', 'seg')
         sex_cmd = ['cex', '%s,%s' % (drz_b, drz_b)]
         sex_cmd += ['-c', os.path.join(self.lr_dir, self.lr_sexfile)]
         # sex_cmd += ['-CATALOG_NAME', '%s.cat' % broot]
         sex_cmd += ['-CATALOG_NAME', '%s.cat' % (os.path.splitext(drz_b)[0])]
         sex_cmd += ['-WEIGHT_IMAGE', '%s,%s' % (wht_b, wht_b)]
         sex_cmd += ['-FLAG_IMAGE', flg_b]
         sex_cmd += ['-CHECKIMAGE_NAME', segimg]
         if len(extra_kwargs):
            for k in extra_kwargs.keys():
               sex_cmd += ['-%s' % k.upper(), str(extra_kwargs[k])]
         print sex_cmd
         subprocess.call(sex_cmd)
      os.chdir(self.lr_dir)

   def resample_psf(self):
      interp_iracpsf.make_iracpsf(self.psffile, 
                     oversample=round(self.lr_scale/self.hr_scale),
                     iracpix=self.lr_scale, radius=150)

   def run_all(self):
      if 'swarp' in self.stages:
         self.swarp()
      elif "blocksum" in self.stages:
         self.blocksum()
         # if blocksum, also shift_crpix
         self.shift_crpix()
      elif "shift_crpix" in self.stages:
         # if running swarp, then no need to run this stage
         self.shift_crpix()
      if 'hst_bgcorr' in self.stages:
         self.hst_bkgd_sub()
      if 'hst_sextractor' in self.stages:
         self.hst_sextractor()
      if 'irac_sextractor' in self.stages:
         self.irac_sextractor()