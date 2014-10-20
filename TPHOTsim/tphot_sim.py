import numpy as np
import fake_galaxies_multires as fgm
from sexsim import sextractor_sim as ss
from sexsim import fake_galaxies
import glob, os, sys
from pygoods import sextractor, Ftable
import yaml
import subprocess
from surfsup import TPHOTpipeline
from pygoods import angsep
from stats import robust, gauss

match_rad = 3.0   # positional match radius (in pixels) for artificial galaxies

class TPHOTSim(ss.SExtractorSim):
   def __init__(self, parfile):
      ss.SExtractorSim.__init__(self, parfile)
      # The first band is the high-res band, and the second band is the 
      # low-res band, and all attributes store the high-res version as the 
      # first element and the low-res attribute as the second element.
      # if self.c.has_key('RA_TARGETS'):
      #    # Only run N galaxies, where N is the number of targets
      #    self.ngal = len(self.c['RA_TARGETS'])
      #    self.ra_targets = self.c['RA_TARGETS']
      #    self.dec_targets = self.c['DEC_TARGETS']
      #    self.target_radius = self.c['TARGET_RADIUS']
      if self.c.has_key('TARGETS'):
         # Only run N galaxies, where N is the number of targets
         self.ngal = len(self.c['TARGETS'])
         self.targetNames = self.c['TARGETS'].keys()
         self.ra_targets = np.array([self.c['TARGETS'][k][1] for k in self.targetNames])
         self.dec_targets = np.array([self.c['TARGETS'][k][2] for k in self.targetNames])
         self.id_targets = np.array([self.c['TARGETS'][k][0] for k in self.targetNames])
         self.target_radius = self.c['TARGET_RADIUS']
      self.tphotparam = self.c['TPHOTPARFILE']
      print "*****************************************"
      print " NUMBER OF SOURCES: %d" % self.ngal
      print "*****************************************"
      fg_args = [self.realimages, self.flagimages, self.bands]
      fg_kwargs = dict(ngal=self.ngal, diskfrac=self.diskfrac, 
                  magfile=self.magfile, rdist=self.rdistfunc, mag0=self.maglow,
                  mag1=self.maghigh, logr0=self.logrmin, logr1=self.logrmax, 
                  lognormal_beta=self.lognormal_beta, 
                  lognormal_mag0=self.lognormal_mag0, 
                  lognormal_peak=self.lognormal_peak,
                  lognormal_sigma=self.lognormal_sigma, 
                  flagmax=self.flagmax,
                  othercols_list=self.othercolnames)
      self.fg = fgm.FakeGalaxiesMultiRes(*fg_args, **fg_kwargs)
      if self.c.has_key('TARGETS'):
         ra, dec = self.fg.get_radec_around(self.ra_targets, self.dec_targets,
                                            radius=self.target_radius)
         x, y = self.fg.get_xy(RA=ra, DEC=dec, mode='hires')
         self.fg.x = x
         self.fg.y = y
      # whether to insert fake sources into the low-res image or not
      self.insert_lores = self.c['INSERT_LORES']  
      self.hires_band = self.bands[0]
      self.lores_band = self.bands[1]
      self.pixscale_hires = self.c['SCALE'][0]
      self.pixscale_lores = self.c['SCALE'][1]

   def insert_fake_sources(self):
      """
      Generate attributes for artificial galaxies, and then insert artificial 
      galaxies into real images.
      """
      # raise NotImplementedError
      if glob.glob('*.list'):
         os.system('rm *.list')
      self.fg.spawn_galaxies()
      if self.c.has_key('TARGETS'):
         ra, dec = self.fg.get_radec_around(self.ra_targets, self.dec_targets,
                                            radius=self.target_radius)
         x, y = self.fg.get_xy(RA=ra, DEC=dec, mode='hires')
         self.fg.x = x
         self.fg.y = y
      if self.insert_lores:
         # Also get the fake sources' (x, y) coordinates in the low-res image
         self.fg.x_lores, self.fg.y_lores = self.fg.get_xy(RA=self.fg.ra, 
                                                     DEC=self.fg.dec, 
                                                     mode='lores')
      # First, add galaxies to the hi-res image
      self.fg.makegals_hires(self.flagimages[self.hires_band])  
      # write input file for mkobjects
      self.makenoiselessimage(self.hires_band, 
                              self.fg.artfiles[self.hires_band], 
                              self.zeropoints[self.hires_band],
                              self.psffiles[self.hires_band])
      self.addsimulated(self.hires_band)
      new_fake_image_hr = os.path.splitext(self.fakeimages[self.hires_band])[0] + '_drz.fits'
      os.rename(self.fakeimages[self.hires_band], new_fake_image_hr)
      self.fakeimages[self.hires_band] = new_fake_image_hr
      # designate a detection image
      self.fake_detect_image = self.fakeimages[self.hires_band]

      # Now insert galaxies into the low-res image
      if self.insert_lores:
         self.fg.makegals_lores()
         self.makenoiselessimage(self.lores_band,
                                 self.fg.artfiles[self.lores_band],
                                 self.zeropoints[self.lores_band],
                                 self.psffiles[self.lores_band])
         self.addsimulated(self.lores_band)
      else:
         self.fakeimages[self.lores_band] = self.realimages[self.lores_band]

   def run_sextractor(self, sex_exe='cex'):
      """
      Run SExtractor in the high-res band.
      Should I also consider running SExtractor in the low-res band?
      """
      assert hasattr(self, 'fg'), "Please generate fake sources first."
      broot = self.root + '_' + self.hires_band
      self.catalogs[self.hires_band] = "%s_%d.cat" % (broot, self.n)
      self.segimages[self.hires_band] = '%s_segnew.fits' % broot         
      self.newcatalogs[self.hires_band] = "%s_%d.newcat" % (broot, self.n)
      n_args = "%s,%s -c %s" % (self.fake_detect_image, self.fakeimages[self.hires_band], self.sexfile)
      n_args = n_args + " -CATALOG_NAME %s" % (self.newcatalogs[self.hires_band])
      n_args = n_args + " -MAG_ZEROPOINT %9.4f" % (self.zeropoints[self.hires_band])
      n_args = n_args + " -GAIN %12.4f" % (self.gains[self.hires_band])
      n_args = n_args + " -FLAG_IMAGE %s" % (self.flagimages[self.hires_band])
      n_args = n_args + " -PARAMETERS_NAME %s" % (self.sextparfile)
      n_args = n_args + " -CHECKIMAGE_TYPE SEGMENTATION"
      n_args = n_args + " -CHECKIMAGE_NAME %s" % (self.segimages[self.hires_band])
      n_args = n_args + " -WEIGHT_TYPE %s,%s" % (self.wht_type, self.wht_type)
      if hasattr(self.c, 'WEIGHT_THRESH'):
         wt = self.c['WEIGHT_THRESH']
         n_args = n_args + " -WEIGHT_THRESH %f,%f" % (wt, wt)
      n_args = n_args + " -WEIGHT_IMAGE %s,%s" % (\
                        self.rmsimages[self.detect_band], self.rmsimages[self.hires_band])
      print "%s %s" % (sex_exe, n_args)
      sys.stdout.flush()
      fpipe = os.popen("%s %s" % (sex_exe, n_args))
      fpipe.close()
      # Identify fake galaxies
      cnew = sextractor(self.newcatalogs[self.hires_band])
      f = open(self.catalogs[self.hires_band], 'wb')
      f.write(cnew._header)
      ncolumns = len(cnew._colnames)
      default_columns = ['X_IN', 'Y_IN', 'MAG_IN', 'RE_IN', 
                        'GTYPE_IN [devauc=%d, expdisk=%d' % (fake_galaxies.devauc, fake_galaxies.disk), 
                        'AXIS_RATIO_IN', 'PA_IN', 'ID_THISRUN']
      for j in range(len(default_columns)):
         ncolumns += 1
         f.write('# %d %s\n' % (ncolumns, default_columns[j]))
      for j in range(len(self.fg.othercolnames)):
         ncolumns += 1
         f.write('# %d %s\n' % ((ncolumns), self.fg.othercolnames[j].upper()))
      if self.magfile:
         f.write('# %d IGAL\n' % (ncolumns+1))
      n_fake_gals = 0
      for i in range(self.ngal):
         dist = np.sqrt((self.fg.x[i]-cnew.x_image)**2 + (self.fg.y[i]-cnew.y_image)**2)
         if dist.min() > match_rad:
            continue
         j = np.argsort(dist)[0]  # index in cnew that matches this fake galaxy
         f.write(' '.join(cnew._colentries[j]))  # write the entry from cnew
         f.write(' %.2f %.2f ' % (self.fg.x[i], self.fg.y[i]))
         f.write(' %.2f %.2f ' % (self.fg.mag[self.hires_band][i], self.fg.re[i]))
         f.write(' %4d  %.2f ' % (self.fg.gtype[i], self.fg.axis_ratio[i]))
         f.write(' %.2f  %4d ' % (self.fg.position_angle[i], i))
         for k in range(len(self.fg.othercolnames)):
            f.write(' %s ' % str(self.fg.othercols[self.fg.othercolnames[k]][i]))
         if self.magfile:
            f.write(' %s ' % str(self.fg.igals[i]))
         f.write('\n')
         n_fake_gals += 1
         self.fg.detected[i] = 1
      # Append non-detected galaxies in the end; substitute all values by -1
      for i in range(self.ngal):
         if self.fg.detected[i] == 0:
            f.write(' '.join([repr(-1)]*len(cnew._colnames)))
            f.write(' %.2f %.2f ' % (self.fg.x[i], self.fg.y[i]))
            f.write(' %.2f %.2f ' % (self.fg.mag[self.hires_band][i], self.fg.re[i]))
            f.write(' %4d  %.2f ' % (self.fg.gtype[i], self.fg.axis_ratio[i]))
            f.write(' %.2f  %4d ' % (self.fg.position_angle[i], i))
            for k in range(len(self.fg.othercolnames)):
               f.write(' %s ' % str(self.fg.othercols[self.fg.othercolnames[k]][i]))
            if self.magfile:
               f.write(' %s ' % str(self.fg.igals[i]))
            f.write('\n')
      f.close()
      print "%d fake galaxies identified." % n_fake_gals
      if n_fake_gals == 0:
         return 0
      os.system("mv %s %s/run%d.cat" % (self.catalogs[self.hires_band], broot, self.n))
      os.system("mv %s %s/run%d.newcat" % (self.newcatalogs[self.hires_band], broot, self.n)) 
      os.system("mv %s %s/glart%d.list" % (self.fg.artfiles[self.hires_band], broot, self.n))                
      os.system("cp %s %s" % (self.psffiles[self.hires_band], broot))

      self._finished_sextractor = True
      return 1

   def write_tphot_param(self, tphotparfile):
      # Only works when self.insert_lores == False
      assert self.insert_lores == False   # should remove this later!!
      broot_hr = self.root + '_' + self.hires_band
      # Pull TPHOT parameters from self.tphotparam?
      tpar = yaml.load(open(self.tphotparam, 'rb'))
      self.tpar = tpar
      f = open(tphotparfile, 'wb')
      f.write('fitdir:     %s\n' % os.getcwd())
      f.write('hiresdir:   %s\n' % tpar['hiresdir'])
      # Also create a symlink within hiresdir that points to the mock image
      curdir = os.getcwd()
      os.chdir(tpar['hiresdir'])
      try:
         os.symlink('%s/%s' % (curdir, self.fakeimages[self.hires_band]), 
                    self.fakeimages[self.hires_band])
      except:
         pass
      try:
         os.symlink('%s/%s' % (curdir, self.segimages[self.hires_band]), 
                    self.segimages[self.hires_band])
      except:
         pass
      if not os.path.exists('%s_%s' % (self.root, self.hires_band)):
         os.symlink('%s/%s_%s' % (curdir, self.root, self.hires_band), 
                    './%s_%s' % (self.root, self.hires_band))
      os.chdir(curdir)
      f.write('loresdir:   %s \n' % tpar['loresdir'])      
      f.write('fit_boxsize:   %.2f \n' % tpar['fit_boxsize'])
      f.write('bkgd_boxsize:  %.2f \n' % tpar['bkgd_boxsize'])
      f.write('growsig:       %.2f \n' % tpar['growsig'])
      f.write('hires_scale:     %.3f \n' % self.pixscale_hires)
      f.write('lores_scale:     %.3f \n' % self.pixscale_lores)
      f.write('fitgauss:        %s \n' % str(tpar['fitgauss']))
      f.write('add_bkgdrms:     %s \n' % str(tpar['add_bkgdrms']))
      f.write('hires_drz:       %s \n' % self.fake_detect_image)
      f.write('hires_seg:       %s \n' % self.segimages[self.hires_band])
      f.write('hires_flag:      %s \n' % self.flagimages[self.hires_band])
      f.write('hires_cat:       %s/run%d.newcat \n' % (broot_hr, self.n))
      f.write('hires_cutoutcat:    run%d_cutout.cat \n' % self.n)
      f.write('hr_fitcat:      %s_run%d_tfit.cat \n' % (self.hires_band.lower(), self.n))
      f.write('hires_band:      %s\n' % self.hires_band.lower())
      ## Again, below only works for self.insert_lores==False!!
      f.write('lores_drz:       %s\n' % self.fakeimages[self.lores_band])
      f.write('lores_seg:       %s\n' % tpar['lores_seg'])  
      # consider re-running SExtractor on low-res image?
      f.write('lores_err:       %s\n' % self.rmsimages[self.lores_band])
      f.write('lores_flg:       %s\n' % self.flagimages[self.lores_band])
      f.write('psffile:         %s\n' % self.psffiles[self.lores_band])
      f.write('lores_band:      %s\n' % self.lores_band.lower())
      f.write('fitpars1:        %s_%s_pass1.param\n' % (self.lores_band.lower(), self.root))
      f.write('fitpars2:        %s_%s_pass2.param\n' % (self.lores_band.lower(), self.root))
      f.write('cutoutdir:       allcut  \n')
      f.write('fitcat:          %s_%s_tphot_pass2.cat\n' % (self.lores_band.lower(), self.root))
      f.write('cluster_name:    %s\n' % self.c['CLUSTER_NAME'])
      f.write('chi2box:         %.1f\n' % tpar['chi2box'])
      f.write('skyrms:          %.4e\n' % tpar['skyrms'])
      # Now figure out which ones are the fake galaxies to run
      c_hr = sextractor('%s_%s/run%d.cat' % (self.root, self.hires_band, self.n))
      # Write target infor
      f.write('targets:\n')
      for i in range(len(c_hr)):
         objInfo = [c_hr.number[i], c_hr.alpha_j2000[i], c_hr.delta_j2000[i]]
         if c_hr.x_image[i] > 0:
            sim_objid = self.n * self.c['ID_OFFSET'] + c_hr.number[i]
            objName = "%s_%d" % (self.root, sim_objid)
            f.write("   %s: %s\n" % (objName, str(objInfo)))
      # write RA
      # f.write('ra:\n')
      # for i in range(len(c_hr)):
      #    if c_hr.x_image[i] > 0:
      #       f.write('   - %.8f \n' % c_hr.alpha_j2000[i])
      # f.write('dec:\n')
      # for i in range(len(c_hr)):
      #    if c_hr.x_image[i] > 0:
      #       f.write('   - %.8f \n' % c_hr.delta_j2000[i])
      # f.write('objnames:\n') 
      # for i in range(len(c_hr)):
      #    if c_hr.x_image[i] > 0:
      #       sim_objid = self.n * self.c['ID_OFFSET'] + c_hr.number[i]
      #       f.write('   - "%s_%d"\n' % (self.root, sim_objid))
      # f.write('objectid:\n')
      # for i in range(len(c_hr)):
      #    if c_hr.x_image[i] > 0:
      #       f.write('   - %d\n' % c_hr.number[i])
      f.close()

   def run_tphot(self, tphotparfile):
      tpipe = TPHOTpipeline.TPHOTpipeline(tphotparfile)
      tpipe.run_all_objects()
      if not self.save:
         # clean up the auxilliary files
         tpipe.clear_aux()
         hires_tfitcat = "%s_run*_tfit.cat" % (self.root)
         os.system("rm %s" % hires_tfitcat)
      return tpipe

   def clear_aux(self, tphotparfile):
      tpipe = TPHOTpipeline.TPHOTpipeline(tphotparfile)
      tpipe.clear_aux()

   def collect_results(self, tpipe):
      """
      Collect the magnitudes from sim sources for a single
      """
      print "*********************************"
      print "Collecting Results..." 
      print "*********************************"
      c_hr = sextractor('%s_%s/run%d.cat' % (self.root, self.hires_band, self.n))
      ngals = len(c_hr)
      # below are the columns from TPHOT catalog
      x_lr = np.zeros(ngals)
      y_lr = np.zeros(ngals)
      Cell = np.zeros(ngals, 'int')
      cx = np.zeros(ngals)
      cy = np.zeros(ngals)
      Rcell = np.zeros(ngals)
      fitqty = np.zeros(ngals)
      fitquerr = np.zeros(ngals)
      sexf = np.zeros(ngals)
      totflux = np.zeros(ngals)
      numfits = np.zeros(ngals, 'int')
      maxflag = np.zeros(ngals, 'int')
      sn = np.zeros(ngals)
      maxcvid = np.zeros(ngals, 'int')
      maxcvratio = np.zeros(ngals)
      mag_lr = np.zeros(ngals)
      magerr_lr = np.zeros(ngals)

      # tpars = yaml.load(open(self.tphotparam))
      
      for i in range(len(c_hr)):
         if c_hr.x_image[i] > 0:
            sim_objid = self.n * self.c['ID_OFFSET'] + c_hr.number[i]
            # print "self.n, ID_OFFSET, numbers[i]", self.n, self.c['ID_OFFSET'], c_hr.number[i]
            # print "sim_objid: %d" % sim_objid
            objname = '%s_%s' % (self.root, sim_objid)
            # print "objname: %s" % objname
            tphot_cat = "%s_%s_tphot_pass2_%s.cat_best_mag" % (self.lores_band.lower(), self.c['CLUSTER_NAME'], objname)
            c_lr = sextractor("%s_%s/" % (self.root, sim_objid) + tphot_cat)
            j = np.arange(len(c_lr))[c_lr.objectid==c_hr.number[i]][0]
            x_lr[i] = c_lr.x[j]
            y_lr[i] = c_lr.y[j]
            Cell[i] = c_lr.cell[j]
            cx[i] = c_lr.cx[j]
            cy[i] = c_lr.cy[j]
            Rcell[i] = c_lr.rcell[j]
            fitqty[i] = c_lr.fitqty[j]
            fitquerr[i] = c_lr.fitquerr[j]
            sexf[i] = c_lr.sexf[j]
            totflux[i] = c_lr.totflux[j]
            numfits[i] = c_lr.numfits[j]
            maxflag[i] = c_lr.maxflag[j]
            sn[i] = c_lr.sn[j]
            maxcvid[i] = c_lr.maxcvid[j]
            maxcvratio[i] = c_lr.maxcvratio[j]
            mag_lr[i] = c_lr.mag[j]
            magerr_lr[i] = c_lr.magerr[j]
      # Now write to an output catalog
      f = open('%s_%s/run%d_comb.cat' % (self.root, self.hires_band, self.n), 'wb')
      f.write(c_hr._header)
      f.write('# %d X_LR\n' % (c_hr._ncolumns + 1))
      f.write('# %d Y_LR\n' % (c_hr._ncolumns + 2))
      f.write('# %d Cell\n' % (c_hr._ncolumns + 3))
      f.write('# %d CX\n' % (c_hr._ncolumns + 4))
      f.write('# %d CY\n' % (c_hr._ncolumns + 5))
      f.write('# %d RCELL\n' % (c_hr._ncolumns + 6))
      f.write('# %d FITQTY\n' % (c_hr._ncolumns + 7))
      f.write('# %d FITQUERR\n' % (c_hr._ncolumns + 8))
      f.write('# %d SEXF \n' % (c_hr._ncolumns + 9))
      f.write('# %d TOTFLUX\n' % (c_hr._ncolumns + 10))
      f.write('# %d NUMFITS\n' % (c_hr._ncolumns + 11))
      f.write('# %d MAXFLAG\n' % (c_hr._ncolumns + 12))
      f.write('# %d SN\n' % (c_hr._ncolumns + 13))
      f.write('# %d MAXCVID\n' % (c_hr._ncolumns + 14)) 
      f.write('# %d MAXCVRATIO \n' % (c_hr._ncolumns + 15))
      f.write('# %d MAG_LR\n' % (c_hr._ncolumns + 16))
      f.write('# %d MAGERR_LR\n' % (c_hr._ncolumns + 17))
      for i in range(len(c_hr)):
         if c_hr.x_image[i] > 0:
            f.write(' '.join(c_hr._colentries[i]))
            f.write(' ')
            f.write('%.2f %.2f ' % (x_lr[i], y_lr[i]))
            f.write('%d %.2f %.2f ' % (Cell[i], cx[i], cy[i]))
            f.write('%.3f %f %f ' % (Rcell[i], fitqty[i], fitquerr[i]))
            f.write('%f %f %d ' % (sexf[i], totflux[i], numfits[i]))
            f.write('%d %f %d ' % (maxflag[i], sn[i], maxcvid[i]))
            f.write('%f %f %f ' % (maxcvratio[i], mag_lr[i], magerr_lr[i]))
            f.write('\n')
      f.close()

   def get_radec(self, objname):
      # return RA, DEC for a target, matching with either the object name 
      # (a string) or object ID (an integer)
      if type(objname) == type('0'):
         try:
            # i = self.c['OBJNAMES'].index(objname)
            i = self.targetNames.index(objname)
         except IndexError:
            print "Object name %s has no match." % objname
            return None
      else:
         try:
            # i = self.c['OBJECTID'].index(objname)
            i = self.id_targets.item(objname)
         except IndexError:
            print "Object ID %d has no match." % objname
            return None
      ra = self.ra_targets[i]
      dec = self.dec_targets[i]
      return (ra, dec)


   def collect_results_all_iter(self, outputfile):
      """
      Merge simulation results from all iterations.
      """
      curdir = os.getcwd()
      os.chdir('%s_%s' % (self.root, self.hires_band))
      cats = glob.glob('run*_comb.cat')
      first = 1
      result_dic = {}
      nobj = 0
      if not len(cats):
         print "No catalogs found."
         return 0
      for cat in cats:
         c = sextractor(cat)
         nobj += len(c)
         if first:
            first = 0
            for d in c._colnames:
               result_dic[d] = c.__getattribute__(d)
         else:
            for d in c._colnames:
               result_dic[d] = np.concatenate([result_dic[d], c.__getattribute__(d)])
      # Write output
      os.chdir(curdir)
      f = open(outputfile, 'wb')
      f.write(c._header)
      for j in range(nobj):
         for i in range(len(c._colnames)):
            f.write('%s ' % str(result_dic[c._colnames[i]][j]))
         f.write('\n')
      f.close()
      # os.chdir(curdir)


class TPHOTSimResults(object):
   """
   Analyze TPHOT simulation results.
   """
   def __init__(self, tphotsim_catalog, ra_col='alpha_j2000', dec_col='delta_j2000', magzero=21.581):
      self.c = sextractor(tphotsim_catalog)
      self.ra_sim = self.c.__getattribute__(ra_col)
      self.dec_sim = self.c.__getattribute__(dec_col)
      self.magzero = magzero

   def calc_nsigma_limits(self, ra, dec, nsigma=3.0, radius=10.0, flux_col='fitqty', print_it=False, full_output=False):
      """
      Calculate the n-sigma magnitude limit within a radius (in arcsec) around
      the supplied sky coordinate (ra, dec). n defaults to 3.
      The simulations are supposed to measure fluxes within the low-res image
      where we did NOT put fake sources... so the spread of the measured flux
      should center around zero and has a dispersion equal to the 1-sigma 
      limiting magnitude. At least in theory...
      """
      angdist = angsep.angsep(ra, dec, self.ra_sim, self.dec_sim)
      pick = (angdist <= (radius / 3600.))
      # fluxerr = self.c.__getattribute__(fluxerr_col)[pick]
      flux = self.c.__getattribute__(flux_col)[pick]
      ### === Use median and MADN to estimate n-sigma flux limit === ###
      median_flux = np.maximum(0., np.median(flux))
      # nsig_fluxerr = np.std(flux) * nsigma  
      # Use MADN for a more robust estimate...
      nsig_fluxerr = robust.MADN(flux) * nsigma
      nsig_flux = median_flux + nsig_fluxerr
      # print "Median flux: ", median_flux
      # print "MADN flux:", nsig_fluxerr / nsigma
      ### ========================================================== ###
      # other ways to estimate disperson? Fit a Gaussian curve?
      ### === Fit Gaussian to the histogram === ###
      # meanflux, sigmaflux = gauss.fitgauss(flux)
      # print "Best-fit Gaussian: mean=%.3e, sigma=%.3e" % (meanflux, sigmaflux)
      # nsig_flux = meanflux + nsigma * sigmaflux
      ### ===================================== ###
      # print "nsig_flux:", nsig_flux
      nsig_mag = self.magzero - 2.5 * np.log10(nsig_flux)
      if print_it:
         print "%.1f-sigma magnitude error from %d simulated sources: %.3f" % (nsigma, pick.sum(), nsig_mag)
      if full_output:
         return nsig_mag, pick.sum(), flux
      else:
         return nsig_mag

