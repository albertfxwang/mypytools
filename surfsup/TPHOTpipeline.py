#!/usr/bin/env python
# Imports a lot of stuff from TFITpipeline.py...

import numpy as np
import TFITpipeline
import os, sys, glob
from pygoods import parseconfig
import archive_tfitrun, tfit_mag
import clean_psf
import subprocess
import yaml
from pygoods import sextractor
import pyfits

class TPHOTpipeline(TFITpipeline.TFITpipeline):
   def __init__(self, paramfile, clean=True, sample='highz'):
      ## Should eventually move the TFITpipeline stuff here if we're making
      ## TFIT obsolete...
      super(TPHOTpipeline, self).__init__(paramfile, clean=clean, sample=sample)
      # This is done in TFITpipeline... I should change how I implement this!
      # if self.c.has_key('targets'):
      #    self.objnames = self.targets.keys()
      #    self.objectid = [self.targets[k][0] for k in self.self.objnames]
      #    self.ra = [self.targets[k][1] for k in self.objnames]
      #    self.dec = [self.targets[k][2] for k in self.objnames]
      # elif self.c.has_key('objectid'):
      # else:
      #    self.objectid = dict(zip(self.objnames, self.c['objectid']))
      os.chdir(self.fitdir)
      print "self.fitdir:", self.fitdir

   def write_tphot_params(self, obj_name):
      """
      Write TPHOT input parameter files.
      """
      os.chdir(self.fitdir)   
      j = np.arange(len(self.ra))[self.objnames==obj_name][0]
      TFITpipeline.cleanfile(self.cutoutdir)
      self.lr_bgcorr = self.lr_root + '_drz_med_%s_bgcorr.fits' % obj_name
      # First for pass 1
      c1 = {}
      c1 = parseconfig('/Users/khuang/Dropbox/codes/tphot/example_pass1.param',
                       paramdict=c1)
      c1['hiresfile'] = os.path.join(self.hiresdir, self.hires_drz)
      c1['hirescat'] = self.hr_fitcat
      c1['hiresseg'] = os.path.join(self.hiresdir, self.hires_seg)
      c1['relscale'] = int(np.round(self.lores_scale/self.hires_scale))
      c1['loresfile'] = os.path.join(self.loresdir, self.lr_bgcorr)
      self.loresfile = c1['loresfile']
      self.lores_err_bkgd = os.path.splitext(self.lores_err)[0] + '_%s_bkgd.fits' % obj_name
      c1['loreserr'] = self.lores_err_bkgd
      c1['kernelfile'] = os.path.join(self.loresdir, self.psffile)
      c1['kernellookup'] = '%s_%s_dancecard_pass1.txt' % (self.lores_band, self.cluster_name)
      c1['cutoutcat'] = '%s_%s_cutout_%s.cat' % (self.lores_band, self.cluster_name, obj_name)
      c1['cutoutdir'] = 'allcut_%s_%s' % (self.lores_band, obj_name)
      # c1['lrfcat'] = '%s_%s_pass1.cat_lrf' % (self.lores_band, self.cluster_name)
      c1['templatedir'] = 'templates_pass1'
      c1['templatecat'] = 'templates_pass1/_templates.cat' 
      c1['fitpars'] = 'tpipe_tphot_pass1.param'
      c1['tphotcat'] = '%s_%s_tphot_pass1.cat' % (self.lores_band, self.cluster_name)
      c1['tphotcell'] = os.path.splitext(c1['tphotcat'])[0] + '.cell'
      c1['tphotcovar'] = os.path.splitext(c1['tphotcat'])[0] + '.covar'
      c1['modelfile'] = '%s_%s_collage_tphot_pass1.fits' % (self.lores_band, self.cluster_name)      
      self.newfitpars1 = os.path.splitext(self.fitpars1)[0]+'_%s.param' % obj_name
      f1 = open(self.newfitpars1, 'wb')
      f1.write('For TPHOT...\n')
      for k in c1.keys():
         if type(c1[k]) == type([]):
            f1.write('%12s %s \n' % (k, ','.join(c1[k])))
         else:
            f1.write('%12s %s \n' % (k, c1[k]))
      f1.close()
      # Then modify the pass 2 parameter file
      c2 = {}
      c2 = parseconfig('/Users/khuang/Dropbox/codes/tphot/example_pass2.param',
                       paramdict=c2)
      c2['hiresfile'] = c1['hiresfile']
      c2['hirescat'] = c1['hirescat']
      c2['hiresseg'] = c1['hiresseg']
      c2['loresfile'] = c1['loresfile']
      c2['loreserr'] = c1['loreserr']
      c2['kernelfile'] = c1['kernelfile']
      c2['kernellookup'] = c1['kernellookup']
      c2['cutoutcat'] = c1['cutoutcat']
      c2['cutoutdir'] = c1['cutoutdir']
      c2['templatedir'] = 'templates_pass2'
      c2['templatecat'] = 'templates_pass1/_templates.cat' 
      c2['tphotpars'] = c1['fitpars'].replace('pass1', 'pass2')
      c2['tphotcat'] = c1['tphotcat'].replace('pass1', 'pass2')
      c2['tphotcell'] = c1['tphotcell'].replace('pass1', 'pass2')
      c2['tphotcovar'] = c1['tphotcovar'].replace('pass1', 'pass2')
      c2['modelfile'] = c1['modelfile'].replace('pass1', 'pass2')
      # if not self.multikern:
      #    c2['multikernels'] = 'false'
      self.newfitpars2 = os.path.splitext(self.fitpars2)[0]+'_%s.param' % obj_name
      f2 = open(self.newfitpars2, 'wb')
      f2.write('For TPHOT...\n')
      for k in c2.keys():
         if type(c2[k]) == type([]):
            f2.write('%12s %s \n' % (k, ','.join(c2[k])))
         else:
            f2.write('%12s %s \n' % (k, c2[k]))
      f2.close()

   def run_tphot(self, obj_name):
      """
      Run TPHOT & clean up and collect the results.
      """
      print "Step 4: Run TPHOT!"
      os.chdir(self.fitdir)
      # Step 4: start running TPHOT
      # os.system('./clean_tphot.sh')
      os.system('rm -r templates_pass?')
      os.chdir(self.loresdir)
      if len(glob.glob('%s_??_??.fits' % self.psfroot)) > 0:
         clean_psf.clean_psf(self.psffile)
      os.chdir(self.fitdir)
      self.newfitpars1 = os.path.splitext(self.fitpars1)[0]+'_%s.param' % obj_name
      self.newfitpars2 = os.path.splitext(self.fitpars2)[0]+'_%s.param' % obj_name
      print "tphot %s" % self.newfitpars1
      subprocess.call(['tphot', self.newfitpars1])
      print "tphot %s" % self.newfitpars2
      subprocess.call(['tphot', self.newfitpars2])
      # After tphot finishes, archive the results
      archive_tfitrun.archive(obj_name)
      # Calculate the magnitudes
      # tphotcatbest = os.path.splitext(self.tphotcat)[0] + '_%s.cat_best' % obj_name
      tphotcatbest = '%s_%s_tphot_pass2_%s.cat_best' % (self.lores_band, self.cluster_name, obj_name)
      tfit_mag.tfit_mag(tphotcatbest)
      os.chdir(self.fitdir)
      os.system('mv %s %s/' % (self.newfitpars1, obj_name))
      os.system('mv %s %s/' % (self.newfitpars2, obj_name))
      os.system('cp %s %s/' % (self.paramfile, obj_name))
      os.system('mv %s %s/' % (self.lr_bgcorr, obj_name))
      self.lr_bgcorr_mask = os.path.splitext(self.lr_bgcorr)[0] + '_mask.fits'
      os.system('mv %s/%s %s/' % (self.loresdir, self.lr_bgcorr_mask, obj_name))
      # os.system('mv %s %s/' % (self.lores_err_bkgd, obj_name))
      os.system('mv chi2mask.fits %s/chi2mask_%s.fits' % (obj_name, obj_name))
      os.chdir(self.hiresdir)
      # os.system('rm %s' % self.hires_seg)
      os.chdir(self.fitdir)
      # os.chdir(obj_name)

   def run_all(self, obj_name):
      """
      Run through all the steps.
      """
      print "*********************************"
      print "Start fitting object %s" % obj_name
      print "*********************************"
      os.chdir(self.fitdir)
      self.hires_catalog(obj_name)
      self.median_bkgd(obj_name)
      self.write_tphot_params(obj_name)
      self.run_tphot(obj_name)
      self.diagnose(obj_name)

   def diagnose(self, obj_name):
      super(TPHOTpipeline, self).diagnose(obj_name)
      os.chdir(self.fitdir)
      if os.path.exists(self.lores_err_bkgd):
         os.system('mv %s %s/' % (self.lores_err_bkgd, obj_name))
      os.system('mv %s %s/' % (self.loresfile, obj_name))
      os.system('mv %s %s/' % (self.lr_bgcorr_mask, obj_name))

   def run_all_objects(self):
      super(TPHOTpipeline, self).run_all_objects()
   
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
         # outcat = os.path.splitext(self.fitcat)[0] + '_%s.cat_best_mag' % obj_name
         outcat = "%s_%s_tphot_pass2_%s.cat_best_mag" % (self.lores_band.lower(), self.cluster_name, obj_name)
         c = sextractor('%s/%s' % (obj_name, outcat))
         mag = c.mag[c.objectid==objid][0]
         mag_err = c.magerr[c.objectid==objid][0]
         flux = c.fitqty[c.objectid==objid][0]
         fluxerr = c.fitquerr[c.objectid==objid][0]
         resid2 = '%s/resid_%s_%s_collage_tphot_pass2_%s.fits' % (obj_name,self.lores_band,self.cluster_name,obj_name)
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

   def get_tphot_mag(self, obj_name):
      if type(obj_name) == type(''):
         objid = self.objectid[self.objnames==obj_name][0]
      else:
         objid = obj_name
         obj_name = self.objnames[self.objectid==objid][0]
      print "SExtractor ID = ", objid
      os.chdir(os.path.join(self.loresdir, 'tphot', obj_name))
      c = sextractor(os.path.splitext(self.fitcat)[0] + '_%s.cat_best_mag' % obj_name)
      mag = c.mag[c.objectid==objid][0]
      magerr = c.magerr[c.objectid==objid][0]
      print "TPHOT magnitudes (w/o aperture correction) for object %s is %.2f +/- %.3f" % (obj_name, mag, magerr)
      resid2 = 'resid_%s_%s_collage_tphot_pass2_%s.fits' % (self.lores_band,self.cluster_name,obj_name)
      hdr = pyfits.getheader(resid2)
      chi2nu = hdr['chi2nu']
      print "Reduced chi2 = %.5f" % chi2nu 
      maxcvratio = c.maxcvratio[c.objectid==objid][0]
      print "Max. covariance ratio: %.3f" % maxcvratio
      os.chdir(os.path.join(self.loresdir, 'tphot'))
      return mag, magerr, chi2nu, maxcvratio

   def get_maxcvratio(self, obj_name):
      if type(obj_name) == type(''):
         objid = self.objectid[self.objnames==obj_name][0]
      else:
         objid = obj_name
         obj_name = self.objnames[self.objectid==objid][0]
      curdir = os.getcwd()
      os.chdir(os.path.join(self.loresdir, 'tphot', obj_name))
      c = sextractor(os.path.splitext(self.fitcat)[0] + '_%s.cat_best_mag' % obj_name)
      os.chdir(curdir)
      return c.maxcvratio[c.objectid==objid][0]

   def get_tphot_mag_all(self):
      mag = np.zeros(len(self.objnames))
      magerr = np.zeros(len(self.objnames))
      for i in range(len(self.objnames)):
         x = self.get_tphot_mag(self.objnames[i])
         mag[i] = x[0]
         magerr[i] = x[1]
      return mag, magerr

   def display_resid(self, obj_name, scales=[-0.01,0.02]):
      os.chdir(os.path.join(self.loresdir, 'tphot', obj_name))
      j = np.arange(len(self.ra))[self.objnames==obj_name][0]
      resid2 = 'resid_%s_%s_collage_tphot_pass2_%s.fits' % (self.lores_band,self.cluster_name,obj_name)
      bgcorr = self.lr_root + '_drz_med_%s_bgcorr.fits' % obj_name
      hires = os.path.join(self.hiresdir, self.hires_drz)
      if os.path.exists(hires):
         cmd = "ds9 %s %s %s -single " % (hires, resid2, bgcorr)
      else:
         cmd = "ds9 %s %s -single " % (resid2, bgcorr)
      cmd += "-regions load all %s.reg " % (obj_name)
      cmd += "-scale limits %f %f " % (scales[0], scales[1])
      cmd += "-match scale "
      cmd += "-pan to %f %f wcs fk5 " % (self.ra[j], self.dec[j])
      cmd += "-zoom to 4 "
      cmd += "-match frame wcs "
      cmd += "&" 
      print cmd
      os.system(cmd)
      ## Also print the residual flux info
      h = pyfits.getheader(resid2)
      print "Reduced chi-square = ", h['chi2nu']
      os.chdir(os.path.join(self.loresdir, 'tphot'))
      

