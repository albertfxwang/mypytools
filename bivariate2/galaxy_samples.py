#!/usr/bin/env python


import numpy as np
import cPickle
from pygoods import Ftable, sextractor
import RLDist
import InterloperDist as ID
import yaml
import pickle_util as pu
import sys, multiprocessing

colnames_default = {'mag': 'mag_out',
   'mag_err': 'mag_out_err',
   're': 're_out_arcsec',
   're_err': 're_out_err_arcsec',
   'n': 'n_out',
   'n_err': 'n_out_err',
   'ar': 'ar_out',
   'ar_err': 'ar_out_err',
   'pa': 'pa_out',
   'pa_err': 'pa_out_err',
   'chi2nu': 'chi2nu',
   'magflag': 'magflag',
   'reflag': 'reflag',
   'nflag': 'nflag',
   'vflag': 'vflag',
   'mag_auto': 'mag_auto',
   'r_hl': 'flux_radius_2_arcsec'}

colnames_old = {'mag': 'magout_gf',
   'mag_err': 'magout_err_gf',\
   're': 'reout_gf',\
   're_err': 'reout_err_gf',\
   'n': 'nout_gf',\
   'n_err': 'nout_err_gf',\
   'ar': 'arout_gf',\
   'ar_err': 'arout_err_gf',\
   'pa': 'paout_gf',\
   'pa_err': 'paout_err_gf',\
   'chi2nu': 'chi2nu_gf',\
   'magflag': 'magflag',\
   'reflag': 'reflag',\
   'nflag': 'nflag',
   'vflag': 'vflag'}

class GalaxySample(object):
   def __init__(self, catalog, field, filtername, z0, maglimits, logrlimits, good_vflags=[0,1,2,3,12], colnames=colnames_default, re_err_lim=0.6, chi2nu_lim=5.0, sexmaglimit=None):
      self.c = Ftable(catalog)
      self.field = field
      self.z0 = z0
      self.colnames = colnames
      self.re_err_lim = re_err_lim
      self.chi2nu_lim = chi2nu_lim
      self.filtername = filtername
      self.maglimits = maglimits
      # if sexmaglimit != None, we will include galaxies fainter than the
      # magnitude limit where GALFIT is still reliable... we include these
      # faint galaxies in order to better constrain the faint-end slope. This
      # is desirable when the GALFIT magnitude limit is not much fainter than
      # L*, e.g., for i-dropouts. In this case, we do not use the size 
      # information at all for these faint galaxies.
      self.sexmaglimit = sexmaglimit
      self.logrlimits = logrlimits
      self.gfsample = np.logical_and(self.check_goodfit(), 
                                     self.check_within())
      self.gfsample = np.logical_and(self.gfsample, self.check_vflag(good_vflags))
      # self.mag, self.logr only include those well-fit by GALFIT and those
      # within maglimits and logrlimits
      self.mag = self.get_col('mag')[self.gfsample]
      self.logr = np.log10(self.get_col('re')[self.gfsample])
      if sexmaglimit == None:
         self.sexsample = np.zeros(len(self.c.d), 'bool')
         self.mag_auto = np.array([])
      else:
         self.sexsample = self.check_within(magcol='mag_auto', recol='r_hl',
                           maglimits=[maglimits[1],sexmaglimit])
         self.mag_auto = self.get_col('mag_auto')[self.sexsample]
      self.data = np.array([self.mag, self.logr]).T  # it is an N by 2 array
      # self._calc_indices = False  # calculate the index of each galaxy on the RL distribution pixel grid

   def check_within(self, magcol='mag', recol='re', maglimits=None, logrlimits=None):
      # Check which objects are within the magnitude and logR limits
      mag_out = self.get_col(magcol)
      re_out = self.get_col(recol)
      logre_out = np.where(re_out > 0, np.log10(re_out), -10.)
      if maglimits == None:
         maglimits = self.maglimits
      if logrlimits == None:
         logrlimits = self.logrlimits
      assert len(maglimits) == 2
      assert len(logrlimits) == 2
      mwithin = np.logical_and(mag_out>=maglimits[0], mag_out<maglimits[1])
      rwithin = np.logical_and(logre_out>=logrlimits[0], 
                               logre_out<logrlimits[1])
      return np.logical_and(mwithin, rwithin)

   def get_vflag(self):
      assert (hasattr(self.c, self.colnames['vflag'])), "Column %s not in catalog %s." % (self.colnames['vflag'],self.c.filename)
      return getattr(self.c, self.colnames['vflag'])

   def get_col(self, attribute):
      colname = self.get_colname(attribute)
      assert (hasattr(self.c, colname)), "Column %s does not exist in the catalog %s." % (colname, self.c.filename)
      return getattr(self.c, colname)

   def get_colname(self, attribute):
      assert (attribute in self.colnames.keys()), "%s not found in column names." % attribute
      return '%s_%s' % (self.filtername, self.colnames[attribute])

   def check_vflag(self, goodflags):
      """
      Check visual inspection flags.
      """
      vflags = self.get_vflag()
      gv = np.in1d(vflags, goodflags)
      print "%d out of %d objects pass visual inspection." % (gv.sum(), len(self.c.d))
      return gv

   def check_goodfit(self):
      # check if GALFIT returned reasonable fits
      # Also enforce magnitude & logR limits
      mag_out = self.get_col('mag')
      recover = (mag_out > 0.)
      re_out = self.get_col('re')
      re_out_err = self.get_col('re_err')
      within = self.check_within()
      good_re = (re_out_err/re_out <= self.re_err_lim)
      good_chi2 = (self.get_col('chi2nu') <= self.chi2nu_lim)
      flags = [self.get_col('magflag'), self.get_col('reflag'), self.get_col('nflag')]
      flagcrit = map(lambda x, y, z: np.any([x,y,z]), *flags)  # if any of the flags are true
      flagcrit = np.logical_not(flagcrit) # only pick the ones with NO flags
      goodfit = reduce(lambda x, y: np.logical_and(x,y), [recover, good_re, good_chi2, flagcrit])
      print 'Number of good fits: %d' % goodfit.sum()
      return goodfit
      # return np.logical_and(goodfit, within)

   def loglikelihood_gfsample(self, M, floor=1.e-50, mag=None, logr=None):
      # Calculate log-likelihood for **a given field**
      # takes M (a RLDist instance) as input
      ### make sure it only uses the GALFIT part of the model distribution ###
      if mag==None: mag = self.mag
      if logr==None: logr = self.logr
      M.value = np.maximum(M.value, floor)
      if self.sexmaglimit != None:
         value = M.window_function(self.maglimits[0], self.maglimits[1],
                                   self.logrlimits[0], self.logrlimits[1])
      else:
         value = M.value.copy()
      # logl_all = np.log(np.array(map(M.__call__, mag, logr)))
      logl_all = np.log(np.array(map(lambda x,y: M.__call__(x,y,value=value),
                      mag, logr)))
      logl = -1. * logl_all.sum()
      return logl

   def loglikelihood_sexsample(self, M, floor=1.e-50):
      # calculate the likelihood component from the faint end (using magnitude
      # number counts only)
      # raise NotImplementedError
      if self.sexmaglimit == None:
         return 0.
      else:
         value = M.window_function(self.maglimits[1], self.sexmaglimit,
                                   self.logrlimits[0], self.logrlimits[1])
         value1d = value.sum(axis=1) * M.dy  # marginalized magnitude dist.
         if value1d.sum() <= floor:
            return 0.
         else:
            sexcounts = np.histogram(self.mag_auto, bins=M.xedges())[0]
            # assert len(value1d) == len(sexcounts)
            # only count the likelihood where the number counts are non-zero
            likelihood = (value1d * sexcounts)[sexcounts>0]
            # print "len(likelihood):", len(likelihood)
            logl = np.log(likelihood) + np.log(M.dx)
            logl_tot = -1. * logl.sum()
            return logl_tot

   def loglikelihood(self, M, floor=1.e-50, mag=None, logr=None):
      logl_gf = self.loglikelihood_gfsample(M, floor=floor, mag=mag, 
                                            logr=logr)
      logl_sex = self.loglikelihood_sexsample(M)  ### to update
      return logl_gf + logl_sex

class LBGSample(GalaxySample):
   # A class that defines the galaxy samples, the limits for each sample, and
   # the transformation kernels for the RL distribution. Parameters here include 
   # everything that is required to calculate the RL distribution, in order to 
   # calculate the likelihood. The method to calculate likelihood is also 
   # defined here.
   def __init__(self, paramfile, filtername, **kwargs):
      # Additional keyword arguments goes to GalaxySample.
      p = yaml.load(open(paramfile))
      self.filtername = filtername
      # The surveyed comoving volume in each field will be factored into the
      # respective dropout-selection kernels, so no need to include them here.
      self.guess = p['GUESS']
      (self.alpha, self.mstar, self.logr0, self.sigma, self.beta) = p['GUESS']
      self.phistar = 0.
      self.fields = p['FIELDS']
      self.z0 = p['Z0']  # nominal redshift of the galaxy sample
      # pixel scales of the image on which the sizes are measured
      self.pixscales = dict(zip(p['FIELDS'], p['PIXSCALES']))  
      self.datafiles = dict(zip(p['FIELDS'], p['DATAFILES']))
      self.maglimits = dict(zip(p['FIELDS'], p['MAGLIMITS']))
      self.logrlimits = dict(zip(p['FIELDS'], p['LOGRLIMITS']))
      if p.has_key('SEXMAGLIMITS'):
         self.sexmaglimits = dict(zip(p['FIELDS'], p['SEXMAGLIMITS']))
         self.sex_kgrids = dict(zip(p['FIELDS'], 
                                map(pu.load_pickle, p['SEX_KGRIDS'])))
      else:
         self.sexmaglimits = dict(zip(p['FIELDS'], [None]*len(self.fields)))
         self.sex_kgrids = dict(zip(p['FIELDS'], [None]*len(self.fields)))
      self.galfit_kgrids = dict(zip(p['FIELDS'], 
                                map(pu.load_pickle, p['GALFIT_KGRIDS'])))
      self.re_err_lim = p['RE_ERR_LIM']
      self.chi2nu_lim = p['CHI2NU_LIM']
      self.technique = p['TECHNIQUE']
      self.magconvert = RLDist.MagConvert(p['MCFILE'])
      self.dropout_kgrids = dict(zip(p['FIELDS'], 
                                 map(pu.load_pickle, p['DROPOUT_KGRIDS'])))
      if p['NPROC_MODEL'] > 0:
         self.nproc_model = p['NPROC_MODEL']
      else:
         self.nproc_model = multiprocessing.cpu_count()
      # interloper fraction parameters
      self.add_interloper = p['ADD_INTERLOPER']
      # self.sdfiles = dict(zip(p['FIELDS'], p['SDFILES']))
      self.size_catalog = p['SIZE_CATALOG']
      self.mag_size_catalog = p['MAG_SIZE_CATALOG']
      self.re_size_catalog = p['RE_SIZE_CATALOG']
      self.intfracfiles = dict(zip(p['FIELDS'], p['INTFRACFILES']))
      self.mag_lolims = dict(zip(p['FIELDS'], p['MAG_LOLIMS']))
      self.goodflags = p['GOOD_VFLAGS']
      if self.add_interloper:
         # Calculate the *multiplicative* interloper fractions in the RL plane
         self.sd = {}
         self.intfrac = {}
         c = Ftable(self.size_catalog)
         mag = getattr(c, self.mag_size_catalog)
         re = getattr(c, self.re_size_catalog)
         for f in self.fields:
            # first calculate the size distribution for all sources
            sd_f=ID.InterloperSizeDist(self.maglimits[f],self.logrlimits[f], 
                                         0.02, 0.02, 0.5, 0.2, f)
            x = sd_f.compute(mag, re)
            self.sd[f] = sd_f
            # then calculate the multiplicative interloper distribution
            intdist=ID.InterloperRLDist(self.maglimits[f],self.logrlimits[f],
                                          0.02, 0.02, f, 
                                          mag_lolim=self.mag_lolims[f])
            intfrac = Ftable(self.intfracfiles[f])
            x = intdist.compute(intfrac, sd_f)
            self.intfrac[f] = intdist

      ### MCMC parameters?
      ### Define a separate class for MCMC stuff...
      self.galaxy_samples = {}
      self.p = p
      # define galaxy samples
      for f in p['FIELDS']:
         self.galaxy_samples[f] = GalaxySample(self.datafiles[f], f, 
                                  self.filtername, self.z0, 
                                  self.maglimits[f], self.logrlimits[f],
                                  re_err_lim=self.re_err_lim,
                                  chi2nu_lim=self.chi2nu_lim,
                                  good_vflags=self.goodflags,
                                  sexmaglimit=self.sexmaglimits[f],
                                  **kwargs)
         print "Total number of galaxies in %s used in likelihood calculation: %d" % (f, len(self.galaxy_samples[f].mag))
      mag, logr = self.combine_samples(flist=self.fields)  # to record attributes self.mag, self.logr
      self.mag = mag
      self.logr = logr
      self.Ntot = len(self.mag)
      # Now builds RLDist factories
      self.RLDist_factories = {}
      factory = RLDist.RLDistributionFactory
      # build RL Distribution factories for each field
      for f in self.fields:
         if self.sexmaglimits[f] == None:
            maglimits = self.maglimits[f]
         else:
            maglimits = [self.maglimits[f][0], self.sexmaglimits[f]]
         self.RLDist_factories[f] = factory(maglimits, 
                                            self.logrlimits[f],
                                            self.magconvert)
      
   def combine_samples(self, flist=None):
      """
      Combine the magnitude and logR from all fiels into one array for mag and 
      one array for logR.
      """
      all_mag = []
      all_logr = []
      if flist == None:
         flist = self.fields
      for f in flist:
         all_mag += [self.galaxy_samples[f].mag]
         all_logr += [self.galaxy_samples[f].logr]
      mag = reduce(lambda x,y: np.concatenate((x,y)), all_mag)
      logr = reduce(lambda x,y: np.concatenate((x,y)), all_logr)
      return mag, logr

   def loglikelihood_tot(self, parameters, floor=1.e-50, verbose=0):
      """
      Calculate log-likelihood given a set of parameters. Returns 
      -log(likelihood).
      Also calculates phi* (as self.phistar).
      """
      # Not doing parallelization here...
      logl_tot = 0.
      model_sum = 0.
      M_list = []
      # first, calculate phistar
      for i in range(len(self.fields)):
         f = self.fields[i]
         M = self.RLDist_factories[f].calculate(parameters, 
               DK=self.dropout_kgrids[f],
               GK=self.galfit_kgrids[f], 
               SK=self.sex_kgrids[f], nproc_model=self.nproc_model,
               gfmaglim=self.maglimits[f][1], semaglim=self.sexmaglimits[f])
         M_list += [M]
         # Each galaxy sample already returns -logl, so no need to multiply by 
         # -1 here again.
         model_sum += M.sum() 
      # calculate phistar, then scale each RL distribution by the *overall*
      # phistar
      self.phistar = float(self.Ntot) / model_sum
      for i in range(len(self.fields)):
         f = self.fields[i]
         # Now calculate the distribution again WITHIN GALFIT limits
         M_list[i] = self.RLDist_factories[f].calculate(parameters, 
               DK=self.dropout_kgrids[f],
               GK=self.galfit_kgrids[f], 
               SK=None, nproc_model=self.nproc_model,
               gfmaglim=self.maglimits[f][1], semaglim=self.maglimits[f][1])
         M_list[i].value = M_list[i].value * self.phistar
         # If add interloper component, DO IT HERE (after normalizing by phistar).
         if self.add_interloper:
            M_list[i].value = M_list[i].value * (1. + self.intfrac[f].value)
         # logl_tot += self.galaxy_samples[f].loglikelihood(M_list[i], 
         #             floor=floor)
         logl_gf = self.galaxy_samples[f].loglikelihood_gfsample(M_list[i],
                     floor=floor)
         logl_tot += logl_gf
         if verbose: print "logl_gf in %s = " % f, logl_gf
         # if we extend the distribution to fainter limits, use SExtractor 
         # kernels
         if self.sexmaglimits[f] != None:
            M_list[i] = self.RLDist_factories[f].apply_transfer_1d(self.sex_kgrids[f], maglim=self.sexmaglimits[f])
            logl_se = self.galaxy_samples[f].loglikelihood_sexsample(M_list[i])
            if verbose: print 'logl_se in %s = ' % f, logl_se
            logl_tot += logl_se
      return logl_tot

   def parameter_string(self, parameters):
      s = ' '.join(map(lambda x: '%.4f'%x, parameters))
      return s

   def mlfunc(self, parameters):
      """
      The function for maximum-likelihood fitting. Just return log-likelihood.
      """
      logl_tot = self.loglikelihood_tot(parameters)
      print "pars:", self.parameter_string(parameters), "logl:", logl_tot
      return logl_tot

   def mlfunc_mcmc(self, parameters):
      """
      For MCMC runs using emcee, also return phistar for incremental recording.
      """
      # Remember to return +1 * log(likelihood)!
      return -1.*self.loglikelihood_tot(parameters), self.phistar

