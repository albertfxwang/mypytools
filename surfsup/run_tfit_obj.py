#!/usr/bin/env python

import numpy as np
import os, sys, subprocess, glob, time
import pyfits
import pywcs
from pyraf import iraf
iraf.stsdas()
import archive_tfitrun
import irac_pipeline
import make_tfitcat
import tfit_mag
import clean_psf
from pygoods import parseconfig  # consider switching to yaml later!!
import yaml
from pygoods import angsep, sextractor

"""
Runs TFIT around a specific object. Estimates the local median background
around the object, then makes necessary masks/catalogs to just fit the objects
surrounding the target. 
This requires a parameter file for the directories of the relevant files as
well as the size of the region to run TFIT on.
At the end of the TFIT run, calculates the magnitude of the target.
Does NOT make PSFs... so user should provide the PSF.
The steps are as follows:
Step 1: makes a cutout catalog in hi-res image that will include the objects
        which will later be included in TFIT.
Step 2: determines local median background level, and subtract it from the 
        low-res science image.
Step 3: re-write TFIT parameter files.
Step 4: run TFIT.
Step 5: collect (archive) this TFIT run.
Step 6: calculates magnitudes and magnitude errors and write the *.cat_best_mag
        catalog.
"""

def cleanfile(filename):
   if os.path.exists(filename):
      if os.path.isdir(filename):
         os.system('rm -r %s' % filename)
      else:
         os.remove(filename)

class TFITpipeline(object):
   def __init__(self, paramfile, clean=True, sample='highz'):
      c = yaml.load(open(paramfile,'rb'))
      self.c = c
      self.paramfile = paramfile
      for k in c.keys():
         setattr(self, k, c[k])
      print self.paramfile
      if not self.hires_drz.endswith('_drz.fits'):
         print "Warning: should rename the high-res science image to the format of *_drz.fits!"
      self.hr_root = os.path.splitext(self.hires_drz)[0][:-4]
      if not self.lores_drz.endswith('_drz.fits'):
         print "Warning: should rename the low-res science image to the format of *_drz.fits!"
      self.lr_root = os.path.splitext(c['lores_drz'])[0][:-4]
      self.objnames = np.array(self.targets.keys())
      self.objectid = np.array([self.targets[k][0] for k in self.objnames])
      self.ra = np.array([self.targets[k][1] for k in self.objnames])
      self.dec = np.array([self.targets[k][2] for k in self.objnames])
      if type(self.bkgd_boxsize) == type(1.0):
         self.bkgd_boxpix = self.bkgd_boxsize / self.lores_scale
      else:
         self.bkgd_boxpix = np.array(self.bkgd_boxsize) / self.lores_scale
      # The size of background box in low-res pixels 
      if type(self.fit_boxsize) == type(1.0):
         self.fit_boxpix = self.fit_boxsize / self.hires_scale
      else:
         self.fit_boxpix = np.array(self.fit_boxsize) / self.hires_scale
      self.hdr_lores = pyfits.getheader(os.path.join(self.loresdir, self.lores_drz))
      self.hdr_hires = pyfits.getheader(os.path.join(self.hiresdir, self.hires_drz))
      # need WCS information in the low-res image to find the object by RA, DEC
      self.hr_wcs = pywcs.WCS(self.hdr_hires)
      self.lr_wcs = pywcs.WCS(self.hdr_lores)
      self.today = time.strftime('%y%m%d')
      self.lores_err_obj = [''] * len(self.targets)
      self.psfroot = os.path.splitext(self.psffile)[0]
      # self.find_matches()
      self.read_positions()
      self.hr_flag_obj = {}
      # dictionaries to store the TFIT results
      self.mag = {}
      self.magerr = {}
      self.sigbkgd = {}
      self.medbkgd = {}
      self.chi2nu = {}
      self.avg_resid = {}
      self.med_resid = {}
      self.sample = sample

   def read_positions(self):
      """
      Read image coordinates for each object.
      """
      os.chdir(self.hiresdir)
      hc = sextractor(self.hires_cat)
      self.x_hires = np.zeros(len(self.ra))
      self.y_hires = np.zeros(len(self.ra))
      self.x_lores = np.zeros(len(self.ra))
      self.y_lores = np.zeros(len(self.ra))
      self.xmin_hr = np.zeros(len(self.ra))
      self.xmax_hr = np.zeros(len(self.ra))
      self.ymin_hr = np.zeros(len(self.ra))
      self.ymax_hr = np.zeros(len(self.ra))
      self.xmin_bkgd = np.zeros(len(self.ra))
      self.xmax_bkgd = np.zeros(len(self.ra))
      self.ymin_bkgd = np.zeros(len(self.ra))
      self.ymax_bkgd = np.zeros(len(self.ra))
      for i in range(len(self.ra)):
         objid = self.objectid[i]
         self.x_hires[i] = hc.x_image[hc.number==objid][0]
         self.y_hires[i] = hc.y_image[hc.number==objid][0]
         lores_xy = self.lr_wcs.wcs_sky2pix([[self.ra[i], self.dec[i]]], 1)[0]
         self.x_lores[i] = lores_xy[0]
         self.y_lores[i] = lores_xy[1]


   def find_matches(self):
      """
      Given the list of RA and DEC, find the ID numbers for the corresponding 
      objects in the high-res catalog.
      """
      os.chdir(self.hiresdir)
      hc = sextractor(self.hires_cat)
      self.objid = np.zeros(len(self.ra), 'int')
      self.match_dist = np.zeros(len(self.ra))  # in arcsec
      self.x_hires = np.zeros(len(self.ra))
      self.y_hires = np.zeros(len(self.ra))
      self.x_lores = np.zeros(len(self.ra))
      self.y_lores = np.zeros(len(self.ra))
      self.xmin_hr = np.zeros(len(self.ra))
      self.xmax_hr = np.zeros(len(self.ra))
      self.ymin_hr = np.zeros(len(self.ra))
      self.ymax_hr = np.zeros(len(self.ra))
      self.xmin_bkgd = np.zeros(len(self.ra))
      self.xmax_bkgd = np.zeros(len(self.ra))
      self.ymin_bkgd = np.zeros(len(self.ra))
      self.ymax_bkgd = np.zeros(len(self.ra))
      for i in range(len(self.ra)):
         angdist = angsep.angsep(self.ra[i], self.dec[i], hc.alpha_j2000, hc.delta_j2000)
         index_min = np.argsort(angdist)[0]
         self.objid[i] = hc.number[index_min]
         self.match_dist[i] = np.min(angdist) * 3600.
         self.x_hires[i] = hc.x_image[index_min]
         self.y_hires[i] = hc.y_image[index_min]
         # Now determine the pixel coordinates in the low-res image
         # The pixel coordinates here are 1-based
         lores_xy = self.lr_wcs.wcs_sky2pix([[self.ra[i], self.dec[i]]], 1)[0]
         self.x_lores[i] = lores_xy[0]
         self.y_lores[i] = lores_xy[1]

   def hires_catalog(self, obj_name, snfloor=0., sub_bkgd=False):
      """
      Make Hi-res TFIT catalog for the cutout stage.
      """
      if obj_name not in self.objnames:
         raise ValueError, "%s not in the list of objects." % obj_name
      print "Step 1: Making hi-res TFIT catalog for cutouts..."
      # variables related to making cutouts
      # Assume that self.hires_drz has the format *_drz.fits
      obj_flag = os.path.splitext(self.hires_flag)[0] + "_%s.fits" % obj_name
      self.hr_flag_obj[obj_name] = obj_flag
      os.chdir(self.hiresdir)      
      cleanfile(obj_flag)
      j = np.arange(len(self.ra))[self.objnames==obj_name][0]
      if not hasattr(self, 'x_hires'):
         print "Have not matched sky coordinates with high-res catalog yet!"
         raise ValueError
      # Below are 1-based pixel indices
      if type(self.fit_boxsize) == type(1.0):
         self.xmin_hr[j] = np.floor(self.x_hires[j] - self.fit_boxpix/2.)
         self.xmax_hr[j] = np.ceil(self.x_hires[j] + self.fit_boxpix/2.)
         self.ymin_hr[j] = np.floor(self.y_hires[j] - self.fit_boxpix/2.)
         self.ymax_hr[j] = np.ceil(self.y_hires[j] + self.fit_boxpix/2.)
      else:
         self.xmin_hr[j], self.xmax_hr[j], self.ymin_hr[j], self.ymax_hr[j] = self.fit_boxsize
      self.xmin_hr[j] = np.maximum(1, int(self.xmin_hr[j]))
      self.xmax_hr[j] = np.minimum(self.hdr_hires['naxis1'], int(self.xmax_hr[j]))
      self.ymin_hr[j] = np.maximum(1, int(self.ymin_hr[j]))
      self.ymax_hr[j] = np.minimum(self.hdr_hires['naxis2'], int(self.ymax_hr[j]))
      print self.xmin_hr[j], self.xmax_hr[j], self.ymin_hr[j], self.ymax_hr[j]
      # Make a flag image for the given target, and convert to 0-based indexing
      flg_hr = pyfits.getdata(self.hires_flag)
      flg_obj_hr = np.ones(flg_hr.shape, 'int32')
      # remember to flip x and y
      flg_obj_hr[self.ymin_hr[j]-1:self.ymax_hr[j],self.xmin_hr[j]-1:self.xmax_hr[j]] = 0
      pyfits.append(obj_flag, flg_obj_hr, self.hdr_hires) 
      # # The old IRAF way
      # iraf.imcalc(c['hires_flag'], self.obj_flag, 
      #          'if (x>=%d)&&(x<%d)&&(y>=%d)&&(y<%d) then 0 else 1' % (c['hr_xmin'], c['hr_xmax'], c['hr_ymin'], c['hr_ymax']),
      #          pixtype='int')
      cleanfile(self.hr_fitcat)
      print os.getcwd()
      make_tfitcat.make_tfitcat(self.hires_cat, flagimage=obj_flag, 
                             segimage=self.hires_seg,
                             outname=os.path.join(self.fitdir, self.hr_fitcat),
                             snfloor=snfloor)
      os.system('rm %s' % obj_flag)
      # set background to zero if sub_bkgd==True; this is the case when one uses
      # a background-subtracted high-res image as input.
      if sub_bkgd:
         hr_fitcat = os.path.join(self.fitdir, self.hr_fitcat)
         c = sextractor(hr_fitcat)
         c.background = 0.
         c.writeto(hr_fitcat, clobber=True)
      os.chdir(self.fitdir)

   def median_bkgd(self, obj_name):
      """
      Determine the local median background in the low-res image, and also add 
      the background uncertainty to the RMS map in quadrature.
      """
      if obj_name not in self.objnames:
         raise ValueError, "%s not in the list of objects." % obj_name
      print "Step 2: determine local median background in low-res images..."   
      os.chdir(self.loresdir)
      j = np.arange(len(self.ra))[self.objnames==obj_name][0]
      lr_bgcorr = self.lr_root + '_drz_med_%s_bgcorr.fits' % obj_name
      self.loresfile = os.path.join(self.loresdir, lr_bgcorr)
      # self.lores_drz = os.path.join(c['loresdir'], c['lores_drz'])
      # self.lores_seg = os.path.join(c['loresdir'], c['lores_seg'])
      # Will include the uncertainty of the local background level in the RMS map
      self.lores_err_obj[j] = os.path.splitext(self.lores_err)[0] + '_%s_bkgd.fits' % obj_name
      # self.loreserr = os.path.join(c['loresdir'], self.loreserr)

      os.chdir(self.fitdir)
      # Now determine the background box
      if type(self.bkgd_boxsize) == type(1.0):
         self.xmin_bkgd[j] = np.floor(self.x_lores[j] - self.bkgd_boxpix)
         self.xmax_bkgd[j] = np.ceil(self.x_lores[j] + self.bkgd_boxpix)
         self.ymin_bkgd[j] = np.floor(self.y_lores[j] - self.bkgd_boxpix)
         self.ymax_bkgd[j] = np.ceil(self.y_lores[j] + self.bkgd_boxpix)
      else:
         self.xmin_bkgd[j], self.xmax_bkgd[j], self.ymin_bkgd[j], self.ymax_bkgd[j] = self.bkgd_boxsize
      self.xmin_bkgd[j] = np.maximum(1, int(self.xmin_bkgd[j]))
      self.xmax_bkgd[j] = np.minimum(self.hdr_lores['naxis1'], int(self.xmax_bkgd[j]))
      self.ymin_bkgd[j] = np.maximum(1, int(self.ymin_bkgd[j]))
      self.ymax_bkgd[j] = np.minimum(self.hdr_lores['naxis2'], int(self.ymax_bkgd[j]))

      subreg_bkgd = [self.xmin_bkgd[j], self.xmax_bkgd[j], self.ymin_bkgd[j], self.ymax_bkgd[j]]
      print "subreg_bkgd:", subreg_bkgd

      #bkgd_boxsize = int(round(c['bkgd_boxsize'] / c['lores_scale']))
      #irac_pipeline.subtract_constant_bkgd(self.lores_drz, self.lores_seg,
      #              self.lr_bgcorr, boxsize=bkgd_boxsize, boxcenter=lr_xypos)
      
      fig = irac_pipeline.subtract_constant_bkgd(os.path.join(self.loresdir, self.lores_drz), 
            os.path.join(self.loresdir, self.lores_seg), 
            os.path.join(self.loresdir, lr_bgcorr), 
            subreg=subreg_bkgd, growsig=self.growsig,
            fitgauss=self.fitgauss)
      fig.axes[0].set_title(fig.axes[0].get_title()+' (ID=%s)' % obj_name)
      if not os.path.isdir(obj_name):
         os.mkdir(obj_name)
      fig.savefig(os.path.join(self.fitdir, '%s/bkgd_hist_%s_%s.png' % (obj_name, obj_name, self.today)))
      hdr_lr_bgcorr = pyfits.getheader(os.path.join(self.loresdir, lr_bgcorr))
      sigbkgd = hdr_lr_bgcorr['SIGBKGD']
      print "Adding local median sky error to the uncertainty map in quadrature..."
      if os.path.exists(self.lores_err_obj[j]):
         os.remove(self.lores_err_obj[j])
      lores_err = pyfits.getdata(os.path.join(self.loresdir, self.lores_err))
      if self.add_bkgdrms:
         lores_err_obj = np.sqrt(lores_err**2 + sigbkgd**2)
      else:
         lores_err_obj = lores_err
      pyfits.append(self.lores_err_obj[j], lores_err_obj, self.hdr_lores)
      self.lr_bgcorr = os.path.join(self.loresdir, lr_bgcorr)
      # iraf.imcalc(os.path.join(c['loresdir'], c['lores_err']), 
      #          os.path.join(c['loresdir'], self.loreserr), 
      #          'sqrt(im1**2 + %f**2)' % self.sigbkgd)

   def write_tfit_params(self, obj_name):
      print "Step 3: Write TFIT parameter files..."
      os.chdir(self.fitdir)   
      j = np.arange(len(self.ra))[self.objnames==obj_name][0]
      cleanfile(self.cutoutdir)
      lr_bgcorr = self.lr_root + '_drz_med_%s_bgcorr.fits' % obj_name
      c1 = {}
      c1 = parseconfig('/Users/khuang/Dropbox/codes/mypytools/tfit_tools/example_tfit_pass1.param',
                       paramdict=c1)
      c1['hiresfile'] = os.path.join(self.hiresdir, self.hires_drz)
      c1['hirescat'] = self.hr_fitcat
      c1['hiresseg'] = os.path.join(self.hiresdir, self.hires_seg)
      c1['relscale'] = int(np.round(self.lores_scale/self.hires_scale))
      c1['loresfile'] = os.path.join(self.loresdir, lr_bgcorr)
      c1['loreserr'] = os.path.splitext(self.lores_err)[0] + '_%s_bkgd.fits' % obj_name
      c1['loresflag'] = os.path.join(self.loresdir, self.lores_flg)
      c1['psffile'] = os.path.join(self.loresdir, self.psffile)
      c1['psflookup'] = '%s_%s_dancecard_pass1.txt' % (self.lores_band, self.cluster_name)
      c1['cutoutcat'] = '%s_%s_cutout_%s.cat' % (self.lores_band, self.cluster_name, obj_name)
      c1['cutoutdir'] = 'allcut_%s_%s' % (self.lores_band, obj_name)
      c1['lrfcat'] = '%s_%s_pass1.cat_lrf' % (self.lores_band, self.cluster_name)
      c1['templatecat'] = '%s_%s_pass1.cat_lrf_culled' % (self.lores_band, self.cluster_name)
      c1['tfitpars'] = 'tpipe_tfit_%s_pass1.param' % self.cluster_name
      c1['tfitcat'] = '%s_%s_tfit_pass1.cat' % (self.lores_band, self.cluster_name)
      c1['tfitcell'] = os.path.splitext(c1['tfitcat'])[0] + '.cell'
      c1['tfitcovar'] = os.path.splitext(c1['tfitcat'])[0] + '.covar'
      c1['modelfile'] = '%s_%s_collage_pass1.fits' % (self.lores_band, self.cluster_name)
      self.newfitpars1 = os.path.splitext(self.fitpars1)[0]+'_%s.param' % obj_name
      f1 = open(self.newfitpars1, 'wb')
      for k in c1.keys():
         if type(c1[k]) == type([]):
            f1.write('%12s %s \n' % (k, ','.join(c1[k])))
         else:
            f1.write('%12s %s \n' % (k, c1[k]))
      f1.close()
      # Then modify the pass 2 parameter file
      c2 = {}
      c2 = parseconfig('/Users/khuang/Dropbox/codes/mypytools/tfit_tools/example_tfit_pass2.param',
                       paramdict=c2)
      c2['hiresfile'] = c1['hiresfile']
      c2['hirescat'] = c1['hirescat']
      c2['hiresseg'] = c1['hiresseg']
      c2['loresfile'] = c1['loresfile']
      c2['loreserr'] = c1['loreserr']
      c2['loresflag'] = c1['loresflag']
      c2['psffile'] = c1['psffile']
      c2['psflookup'] = c1['psflookup']
      c2['cutoutcat'] = c1['cutoutcat']
      c2['cutoutdir'] = c1['cutoutdir']
      c2['lrfcat'] = c1['lrfcat']
      c2['templatecat'] = '%s_%s_pass2.cat_lrf_culled' % (self.lores_band, self.cluster_name)
      c2['tfitpars'] = c1['tfitpars'].replace('pass1', 'pass2')
      c2['tfitcat'] = c1['tfitcat'].replace('pass1', 'pass2')
      c2['tfitcell'] = c1['tfitcell'].replace('pass1', 'pass2')
      c2['tfitcovar'] = c1['tfitcovar'].replace('pass1', 'pass2')
      c2['modelfile'] = c1['modelfile'].replace('pass1', 'pass2')
      self.newfitpars2 = os.path.splitext(self.fitpars2)[0]+'_%s.param' % obj_name
      f2 = open(self.newfitpars2, 'wb')
      for k in  c2.keys():
         if type(c2[k]) == type([]):
            f2.write('%12s %s \n' % (k, ','.join(c2[k])))
         else:
            f2.write('%12s %s \n' % (k, c2[k]))
      f2.close()

   def run_tfit(self, obj_name):
      """
      Run TFIT & clean up and collect the results.
      """
      print "Step 4: Run TFIT!"
      os.chdir(self.fitdir)
      # Step 4: start running TFIT
      os.system('./clean_tfit.sh')
      os.system('rm -r templates_pass?')
      os.chdir(self.loresdir)
      if len(glob.glob('%s_??_??.fits' % self.psfroot)) > 0:
         clean_psf.clean_psf(self.psffile)
      os.chdir(self.fitdir)
      self.newfitpars1 = os.path.splitext(self.fitpars1)[0]+'_%s.param' % obj_name
      self.newfitpars2 = os.path.splitext(self.fitpars2)[0]+'_%s.param' % obj_name
      print "tfit %s" % self.newfitpars1
      subprocess.call(['tfit', self.newfitpars1])
      print "tfit %s" % self.newfitpars2
      subprocess.call(['tfit', self.newfitpars2])
      # After TFIT finishes, archive the results
      archive_tfitrun.archive(obj_name)
      # Calculate the magnitudes
      # tfitcatbest = os.path.splitext(self.tfitcat)[0] + '_%s.cat_best' % obj_name
      tfitcatbest = '%s_%s_tfit_pass2_%s.cat_best' % (self.lores_band, self.cluster_name, obj_name)
      tfit_mag.tfit_mag(tfitcatbest)
      os.chdir(self.fitdir)
      os.system('mv %s %s/' % (self.newfitpars1, obj_name))
      os.system('mv %s %s/' % (self.newfitpars2, obj_name))
      os.system('cp %s %s/' % (self.paramfile, obj_name))
      os.chdir(obj_name)

   def diagnose(self, obj_name):
      """
      Run diagnostics on how TFIT runs.
      """
      os.chdir(self.fitdir)
      print "****************************"
      print "Diagnosing %s " % obj_name
      print "****************************"
      # Calculate chi-square around the target
      self.newfitpars2 = os.path.splitext(self.fitpars2)[0]+'_%s.param' % obj_name
      self.newfitpars2 = '%s/%s' % (obj_name, self.newfitpars2)
      fitpars2 = {}
      fitpars2 = parseconfig(self.newfitpars2, paramdict=fitpars2)
      newcollage = os.path.splitext(fitpars2['modelfile'])[0] + '_%s.fits' % obj_name
      newcollage = '%s/%s' % (obj_name, newcollage)
      collage = pyfits.getdata(newcollage)
      loresfile_data = pyfits.getdata(fitpars2['loresfile'])
      loreserr_data = pyfits.getdata(fitpars2['loreserr'])
      j = np.arange(len(self.ra))[self.objnames==obj_name][0]
      xypos = [self.x_lores[j], self.y_lores[j]]
      chi2box = [xypos[0]-self.chi2box/2-1, xypos[0]+self.chi2box/2-1, 
                 xypos[1]-self.chi2box/2-1, xypos[1]+self.chi2box/2-1]
      # print chi2box
      chi2mask = np.ones(collage.shape, 'int')
      # print np.shape(chi2mask)
      chi2mask[chi2box[2]:chi2box[3]+1,chi2box[0]:chi2box[1]+1] = 0
      collage_masked = np.ma.array(collage, mask=chi2mask)
      loresfile_masked = np.ma.array(loresfile_data, mask=chi2mask)
      chi2 = np.ma.sum(((collage_masked - loresfile_data) / loreserr_data)**2)
      chi2nu = chi2 / np.sum(chi2mask == 0)
      residfile = os.path.splitext(fitpars2['modelfile'])[0] + '_%s.fits' % obj_name
      residfile = '%s/resid_%s' % (obj_name, residfile)
      resid = np.ma.array(pyfits.getdata(residfile), mask=chi2mask)
      mean_resid = np.ma.mean(resid)
      median_resid = np.ma.median(resid)

      hdu_chi2mask = pyfits.PrimaryHDU(chi2mask)
      cleanfile('chi2mask.fits')
      hdu_chi2mask.writeto('chi2mask.fits')
      f = open('%s/%s_diags_tfit.txt' % (obj_name, obj_name), 'wb')
      print >> f, "Total chi-square over %d pixels:" % np.sum(chi2mask==0), chi2
      print >> f, "Reduced chi-square: ", chi2nu
      print >> f, "Mean residual flux within the box: %f (%.2f%% of sky RMS)" % (mean_resid, 100.*mean_resid/self.skyrms)
      print >> f, "Median residual flux within the box: %f (%.2f%% of sky RMS" % (median_resid, 100.*median_resid/self.skyrms)
      f.close()

      h_resid = pyfits.open(residfile, mode='update')
      h_resid[0].header['LRFILE'] = os.path.split(fitpars2['loresfile'])[-1]
      h_resid[0].header['LRERR'] = os.path.split(fitpars2['loreserr'])[-1]
      h_resid[0].header['PSFFILE'] = self.psffile
      h_resid[0].header['XMIN_HR'] = self.xmin_hr[j]
      h_resid[0].header['XMAX_HR'] = self.xmax_hr[j]
      h_resid[0].header['YMIN_HR'] = self.ymin_hr[j]
      h_resid[0].header['YMAX_HR'] = self.ymax_hr[j]
      h_resid[0].header['XMINBKGD'] = self.xmin_bkgd[j]
      h_resid[0].header['XMAXBKGD'] = self.xmax_bkgd[j]
      h_resid[0].header['YMINBKGD'] = self.ymin_bkgd[j]
      h_resid[0].header['YMAXBKGD'] = self.ymax_bkgd[j]
      h_resid[0].header['CHI2DOF'] = np.sum(chi2mask==0)
      h_resid[0].header['CHI2TOT'] = chi2
      h_resid[0].header['CHI2NU'] = chi2nu
      h_resid[0].header['AVGRESID'] = mean_resid
      h_resid[0].header['MEDRESID'] = median_resid
      h_resid.flush()
      h_resid.close()
   
      # Make cutouts of the input and residual images in the low-res frame.
      lr_bgcorr = self.lr_root + '_drz_med_%s_bgcorr.fits' % obj_name
      lr_bgcorr = os.path.join(self.loresdir, lr_bgcorr)
      hdr = pyfits.getheader(lr_bgcorr)
      bgcorr = pyfits.getdata(lr_bgcorr)
      resid = pyfits.getdata(residfile)
      # if type(self.fit_boxsize) == type(1.0):
      #    fit_boxpix_lr = self.fit_boxsize / self.lores_scale
      #    xmin_lr = int(np.floor(self.x_lores[j]-fit_boxpix_lr/2.))
      #    xmax_lr = int(np.ceil(self.x_lores[j]+fit_boxpix_lr/2.))
      #    ymin_lr = int(np.floor(self.y_lores[j]-fit_boxpix_lr/2.))
      #    ymax_lr = int(np.floor(self.y_lores[j]+fit_boxpix_lr/2.))
      # else:
      #    xmin_lr, ymin_lr = 
      ra_min, dec_min = self.hr_wcs.wcs_pix2sky([[self.xmin_hr[j], self.ymin_hr[j]]], 1)[0]
      xmin_lr, ymin_lr = self.lr_wcs.wcs_sky2pix([[ra_min, dec_min]], 1)[0]
      ra_max, dec_max = self.hr_wcs.wcs_pix2sky([[self.xmax_hr[j], self.ymax_hr[j]]], 1)[0]
      xmax_lr, ymax_lr = self.lr_wcs.wcs_sky2pix([[ra_max, dec_max]], 1)[0]
      hdr['crpix1'] = hdr['crpix1'] - (xmin_lr-1)
      hdr['crpix2'] = hdr['crpix2'] - (ymin_lr-1)
      cleanfile('%s/%s_lr_cut_%s.fits' % (obj_name, self.lores_band, obj_name))
      cleanfile('%s/%s_resid_%s.fits' % (obj_name, self.lores_band, obj_name))
      pyfits.append('%s/%s_lr_cut_%s.fits' % (obj_name, self.lores_band, obj_name),
                    bgcorr[ymin_lr-1:ymax_lr,xmin_lr-1:xmax_lr], hdr)
      pyfits.append('%s/%s_resid_%s.fits' % (obj_name, self.lores_band, obj_name),
                    resid[ymin_lr-1:ymax_lr,xmin_lr-1:xmax_lr], hdr)
      # write a region file for this object
      f = open('%s/%s.reg' % (obj_name, obj_name), 'wb')
      f.write('fk5; circle(%f, %f, 1.0") #text={%s}\n' % (self.ra[j], self.dec[j], obj_name))
      f.close()
      # Now collect results
      x = glob.glob('%s/*.cat_best_mag' % obj_name)[0]
      c = sextractor(x)
      self.mag[obj_name] = c.mag[c.objectid==self.objectid[j]][0]
      self.magerr[obj_name] = c.magerr[c.objectid==self.objectid[j]][0]
      self.sigbkgd[obj_name] = h_resid[0].header['sigbkgd']
      self.medbkgd[obj_name] = h_resid[0].header['medbkgd']
      self.chi2nu[obj_name] = chi2nu
      self.avg_resid[obj_name] = mean_resid
      self.med_resid[obj_name] = median_resid
      
   def diagnose_all(self):
      for obj in self.objnames:
         self.diagnose(obj)

   def run_all(self, obj_name):
      """
      Run through all the steps.
      """
      self.hires_catalog(obj_name)
      self.median_bkgd(obj_name)
      self.write_tfit_params(obj_name)
      self.run_tfit(obj_name)
      self.diagnose(obj_name)

   def run_all_objects(self):
      """ 
      Run through all objects.
      """
      for obj_name in self.objnames:
         self.run_all(obj_name)

   def collect_results(self):
      # First get the MAG_ISO from the high-res catalog
      ch = sextractor(os.path.join(self.hiresdir,self.hires_cat))
      mag_iso = {}
      magerr_iso = {}
      f = open('%s_%s_tfit_%s.cat' % (self.lores_band, self.cluster_name, self.sample), 'wb')
      f.write('# 1 OBJNAME\n')
      f.write('# 2 RA\n')
      f.write('# 3 DEC\n')
      f.write('# 4 %s_MAG\n' % self.lores_band)
      f.write('# 5 %s_MAG_ERR\n' % self.lores_band)
      f.write('# 6 %s_MED_BKGD\n' % self.lores_band)
      f.write('# 7 %s_SIG_BKGD\n' % self.lores_band)
      f.write('# 8 %s_AVG_RESID\n' % self.lores_band)
      f.write('# 9 %s_MED_RESID\n' % self.lores_band)
      f.write('# 10 %s_CHI2NU\n' % self.lores_band)
      f.write('# 11 %s_MAG_ISO\n' % self.hires_band)
      f.write('# 12 %s_MAGERR_ISO\n' % self.hires_band)
      f.write('# 13 %s_FLUX\n' % self.lores_band)
      f.write('# 14 %s_FLUXERR\n' % self.lores_band)
      for obj_name in self.objnames:
         j = np.arange(len(self.ra))[self.objnames==obj_name][0]
         objid = self.objectid[j]
         mag_iso[obj_name] = ch.mag_iso[ch.number==objid][0]
         magerr_iso[obj_name] = ch.magerr_iso[ch.number==objid][0]
         outcat = os.path.splitext(self.fitcat)[0] + '_%s.cat_best_mag' % obj_name
         c = sextractor('%s/%s' % (obj_name, outcat))
         mag = c.mag[c.objectid==objid][0]
         mag_err = c.magerr[c.objectid==objid][0]
         flux = c.fitqty[c.objectid==objid][0]
         fluxerr = c.fitquerr[c.objectid==objid][0]
         resid2 = '%s/resid_%s_%s_collage_pass2_%s.fits' % (obj_name,self.lores_band,self.cluster_name,obj_name)
         hdr = pyfits.getheader(resid2)
         medbkgd = hdr['medbkgd']
         sigbkgd = hdr['sigbkgd']
         avgresid = hdr['avgresid']
         medresid = hdr['medresid']
         chi2nu = hdr['chi2nu']
         f.write('%-15s %f %f ' % (obj_name, self.ra[j], self.dec[j]))
         f.write('%.3f %.3f ' % (mag, mag_err))
         f.write('%f %f ' % (medbkgd, sigbkgd))
         f.write('%f %f ' % (avgresid, medresid))
         f.write('%f ' % chi2nu)
         f.write('%f %f ' % (mag_iso[obj_name], magerr_iso[obj_name]))
         f.write('%f %f ' % (flux, fluxerr))
         f.write('\n')
      f.close()

   def clear_all(self):
      self.clear_aux()
      for obj_name in self.objnames:
         os.system('rm -r %s' % obj_name)

   def clear_aux(self):
      os.system('rm -r templates_pass*')
      os.system('rm -r multikern')
      os.system('rm tpipe*')
      os.system('rm FluxOrderedCat.txt')
      os.system('rm *_pass?.cat')
      os.system('rm *.covar')
      os.system('rm *dancecard*')
      for obj_name in self.objnames:
         os.system('rm -r covardir*')
         os.system('rm -r allcut_*_%s' % obj_name)
         os.system('rm -r *_cutout_%s.cat*' % obj_name)
         os.system('rm -r *_unc_%s_bkgd.fits' % obj_name)
         os.system('rm %s/*_flg_%s.fits' % (self.hiresdir, obj_name))
         os.system('rm -r covardir*')

class run_tfit_obj(TFITpipeline):
   # An alias for the older version of the class name
   def __init__(self, paramfile, clean=True, sample='highz'):
      super(run_tfit_obj, self).__init__(paramfile, clean=clean, sample=sample)


def merge_catalogs(catalogs, outname):
   catlist = []
   colnames = ['objname', 'ra', 'dec']
   columns = []
   for cat in catalogs:
      c = sextractor(cat)
      if not len(columns):
         columns += [c.objname]
         columns += [c.ra]
         columns += [c.dec]
      for col in c._colnames:
         if col not in colnames:
            colnames += [col]
            columns += [c.__getattribute__(col)]
   print colnames
   print np.shape(columns)
   f = open(outname, 'wb')
   for i in range(len(colnames)):
      f.write('# %d %s\n' % ((i+1), colnames[i].upper()))
   for j in range(len(c.ra)):
      for i in range(len(colnames)):
         f.write('%s ' % columns[i][j])
      f.write('\n')
   f.close()

if __name__ == "__main__":
   paramfile = sys.argv[1]
   r = run_tfit_obj(paramfile)
   r.run_all()
