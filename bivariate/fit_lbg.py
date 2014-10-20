#!/usr/bin/env python

from numpy import *
#import mlutil
from pygoods import *
import os, sys
import bivRL as bl
import bivariate_fit as bf
import time
import copy
import matplotlib.pyplot as plt
#import fig_library as fl
import interlopers
import scipy
from scipy import stats
import cPickle
import zdist
import emcee
import operator
import interlopers
import yaml

__version__ = "2.0 22 May 2013"
__author__ = "Kuang-Han Huang; JHU/STScI"

today = time.strftime('%y%m%d')

"""
Revision History:
v2.0  -- revised to fit an arbitrary number of data sets
"""

def readkgrid(kgridfile):
   f = open(kgridfile, 'rb')
   kgrid = cPickle.load(f)
   f.close()
   return kgrid

def replaceNone(a):
   """Works only for 1-D arrays/lists"""
   for i in range(len(a)):
      if a[i]=="None":
         a[i] = None
   return a

def between(arr, limits):
   x = (operator.ge(arr,limits[0]) & operator.lt(arr,limits[1]))
   return x



## Needs revision!! (130522)
def quickmlfunc(par, dlimits1, dlimits2, kgrid1, kgrid2, zdgrid1, zdgrid2, mc,
   drop='b', verbose=0):
   if drop == 'b': 
      zmean = 4.0
      lbgcat1 = 'bdrops_gf_v2.cat'
      lbgcat2 = 'bdrops_udf_gf_v2.cat'
      chisqnulim1 = 0.4
   elif drop == 'v': 
      zmean = 5.0
      lbgcat1 = 'vdrops_gf_v2.cat'
      lbgcat2 = 'vdrops_udf_gf_v2.cat'
      chisqnulim1 = 0.5
   mag1, re1, crit1 = cleandata(lbgcat1, reerr_ratiolim=0.6, chisqnulim=chisqnulim1,\
      cullflags=[0,1,3,12],limits=dlimits1, drop=drop) 
   # GOODS data
   data1 = array([mag1, log10(re1)])
   mag2, re2, crit2 = cleandata(lbgcat2, reerr_ratiolim=0.6, chisqnulim=5.0,\
      cullflags=[0,1,3,12], limits=dlimits2, drop=drop)
   data2 = array([mag2, log10(re2)])
   data2 = array([[],[]])
   #print "shape(data1)", shape(data1)
   #print "shape(data2)", shape(data2)
   #data = concatenate((data1, data2), axis=1)
   # UDF data
   logl = bf.mlfunc(par, data1, data2, dlimits1, dlimits2, bl.pixdx,\
      kgrid1, kgrid2, 1.0, 1.0, -21.0, zdgrid1, zdgrid2, mc, \
      1, mc(zmean), verbose, -1, 'phistar', drop)
   return logl

class FitBivariateRL(bf.biv_fit): 
   # switching to the object-oriented approach
   # a base class for fitting size-luminosity distributions
   # Does NOT take care of data initialization (e.g., cleaning of interlopers)
   def __init__(self, parfile):
      # reads in the parameter file that contains everything
      #meanz_dic = {'b':4.0, 'v':5.0}
      self.parfile = parfile
      #c = parseconfig(parfile)
      c = yaml.load(open(parfile, 'rb'))
      for k in c.keys():
         print k, '=', c[k]
      for k in c.keys():
         setattr(self, k, c[k])
         # set all entries in the parameter file as attributes
      self.parameters = array([self.ALPHA, self.MSTAR, self.LOGR0, self.SIGMA, self.BETA])
      # IMPORTANT:
      # WHENEVER A FITTING IS PERFORMED, ONE NEEDS TO INITIALIZE
      # THE FOLLOWING ATTRIBUTES FOR THE FITTING TO WORK:
      # self.datasets
      # self.tfkgrid
      # self.zdgrid
      # self.limits
      # These all should be dictionaries and cannot be easily included
      # in the parameter files. However, here I use initiliaze them to 
      # null values.
      self.datasets = {}
      self.tfkgrid = {}
      self.zdgrid = {}
      self.limits = {}
      for f in self.FIELDS:
         self.datasets[f] = None
         self.tfkgrid[f] = None
         self.zdgrid[f] = None
         self.limits[f] = None
      self.nfields = len(self.FIELDS)   
      self.pixdx = self.PIXDX
      self.drop = self.DROP
      if hasattr(self, 'ALPHA'):
         self.guess = array([self.ALPHA, self.MSTAR, self.LOGR0, self.SIGMA,
                            self.BETA])
      self.M0 = self.M0
      self.z_mean = self.Z_MEAN
      if hasattr(self, 'MCFILE'):
         # Same MCFILE for all fields
         mc0 = bl.mconvert(self.MCFILE)
         self.mc = {}
         for f in self.FIELDS:
            self.mc[f] = mc0

      self.xout = array([])
      # initialize arguments for bf.mlfunc
      self.fft = 1  # use FFT for convolution
      self.norm = -1.  # Do not normalize the model b/c we'll use phistar
      self.technique = self.TECHNIQUE
      self.phistar = []
      self.p_history = []
      # models
      self.sd = {}
      self.intfrac = {}
      self.mag_lolims = {}
      # if self.ADD_INTERLOPER.lower() in ['yes','true']:
      self.add_interloper = self.ADD_INTERLOPER
      #self.sdfiles = self.SDFILES
      #self.intfrac_files = self.INTFRAC_FILES
      if self.add_interloper == False:
         for f in self.FIELDS:
            self.sd[f] = None
            self.intfrac[f] = None
            self.mag_lolims[f] = 21.0
      self.nproc_model = self.NPROC_MODEL
      # initialize the MCMC-related parameters
      self.nwalkers = self.NWALKERS
      self.nthreads_mcmc = self.NTHREADS_MCMC
      self.nburn_mcmc	 = self.NBURN_MCMC
      self.niter_mcmc = self.NITER_MCMC
      #print "Total points used in fitting:", (len(mag1)+len(mag2))
      ####
      self.sampler = None
      self.verbose = 1

   def calc_phistar(self, params):
      self.mlfunc_RL(params, *self.mlfunc_args)
      return self.phistar[-1]

   def mlfit(self):
      # designate self.mlfunc
      self.mlfunc = self.mlfunc_RL
      print "Start at:", time.ctime()
      t1 = time.time()
      # entering fitting code
      print "Start optimization"
      xout = self.bivariate_fit(self.mlfunc_args, self.guess, self.drop,
                              self.boundary, technique=self.technique)
      t2 = time.time()
      xout[3] = abs(xout[3]); xout[4] = abs(xout[4])
      print "After fitting:", time.ctime()
      print "Duration:", (t2-t1), "seconds"
      sys.stdout.flush()
      print "Best-fit parameters:"
      print xout
      self.xout = xout
      self.phistar = array(self.phistar)
      self.p_history = array(self.p_history)


   def build_sampler(self, nwalkers, ndim=5):
      # use one thread for now... maybe figure out later how to use
      # multiprocessing
      # use emcee to perform MCMC sampling of the posterior probability distribution
      self.nwalkers = nwalkers
      self.sampler = emcee.EnsembleSampler(nwalkers, ndim, self.mlfunc_RL_mcmc, 
                                   args=self.mlfunc_args, threads=1)
      self.has_walkers = False

   def create_walkers(self, par0, par0_scale=[0.1,0.5,0.1,0.1,0.1]):
      alpha0 = random.uniform(par0[0]-par0_scale[0], 
                                 par0[0]+par0_scale[0], 
                                 size=self.nwalkers)
      mstar0 = random.uniform(par0[1]-par0_scale[1], 
                                 par0[1]+par0_scale[1], 
                                 size=self.nwalkers)
      logr00 = random.uniform(par0[2]-par0_scale[2], 
                              par0[2]+par0_scale[2], 
                              size=self.nwalkers)
      sigma0 = random.uniform(par0[3]-par0_scale[3], 
                                 par0[3]+par0_scale[3], 
                                 size=self.nwalkers)
      beta0 = random.uniform(par0[4]-par0_scale[4], 
                                par0[4]+par0_scale[4], 
                                size=self.nwalkers)
      par0 = concatenate([alpha0,mstar0,logr00,sigma0,beta0])
      par0 = par0.reshape(5,self.nwalkers).swapaxes(0,1)
      self.par0_walkers = par0
      self.has_walkers = True

   def mcmc_sample(self, par0, par0_scale=[0.1,0.5,0.1,0.1,0.1], 
                   iterations=1, progress_file=None, verbose=0):
      # sample the posterior PDF using MCMC
      #if self.nthreads_mcmc > 1:
      #   self.nproc_model = 1
      #   self.set_mlfunc()
      self.verbose = verbose
      # Modified mlfunc_RL so that if self.sampler!=None, mlfunc_RL will
      # return two numbers, loglikehood and phistar. The phistar will be
      # stored in self.blobs, as shape (niterations, nwalkers).
      if self.sampler == None:
         print "MCMC sampler not built yet."
         return 0
      else:
         # first randomly generates the initial positions of the 
         # walkers around the initial guess
         if progress_file == None:
            progress_file = "mcmc_run_%s.fits" % today
         # incrementally save MCMC progress, following the recipe on the EMCEE 
         # website.
         if not os.path.exists(progress_file):
            f = open(progress_file, 'w')
            f.close()
         if self.has_walkers == False:
            self.create_walkers(par0, par0_scale)
            # pos, prob, state, blobs = self.sampler.run_mcmc(self.par0_walkers, 
            #                                              iterations)
            for result in self.sampler.sample(self.par0_walkers, 
                                              iterations=iterations, 
                                              storechain=False):
               position = result[0]
               blobs = result[3]
               f = open(progress_file, 'a')
               accept_ratio = self.sampler.acceptance_fraction
               for k in range(position.shape[0]):
                  # print "k=", k
                  # f.write('%d %s\n' % (k, " ".join(position[k])))
                  f.write('%d %s ' % (k, " ".join(map(lambda x: "%f"%x, position[k]))))
                  f.write('%.8f %.3f \n' % (blobs[k], accept_ratio[k]))
               f.close()
            pos, prob, state, blobs = result
         else:
            pos, prob, state, blobs = self.sampler.run_mcmc(self.mcmc_pos, 
                                                            iterations)
         self.mcmc_pos = pos
         self.mcmc_prob = prob
         self.mcmc_state = state
   
   def write_mcmc(self, catalog):
      f = open(catalog, 'wb')
      f.write('# 1 ALPHA\n')
      f.write('# 2 MSTAR\n')
      f.write('# 3 LOGR0\n')
      f.write('# 4 SIGMA\n')
      f.write('# 5 BETA\n')
      f.write('# 6 PHISTAR\n')
      flatchain = self.sampler.flatchain
      blobs = array(self.sampler.blobs).ravel()
      for i in range(len(blobs)):
         chain = flatchain[i]
         f.write('%.4f %.4f %.4f %.4f %.4f ' % (chain[0],chain[1],chain[2],
                 chain[3],chain[4]))
         f.write('%.4f ' % blobs[i])
         f.write('\n')
      f.close()

   def pickle(self, pname):
      # save the current state
      f = open(pname, 'w')
      cPickle.dump(self, f, 2)
      f.close()

class FitBivariateRM(FitBivariateRL):
   """
   A super class of FitBivariateRL that performs fitting of the size-stellar
   mass distribution. Just override the mlfit method.
   """
   #def __init__(self, parfile):
   #   FitBivariateRL.__init__(self, parfile)
   #   self.phistar_m = []
   #   self.logM0 = sefl.LOGM0

   def mlfit(self):
      # specify self.mlfunc
      self.mlfunc = self.mlfunc_RM
      self.guess = array([self.ALPHA_M, self.LOGMSTAR_M, self.LOGR0, 
                         self.SIGMA, self.BETA_M])
      print "self.mlfunc", self.mlfunc
      print "Start at:", time.ctime()
      t1 = time.time()
      # entering fitting code
      print "Start optimization"
      # need work here
      xout = self.bivariate_fit(self.mlfunc_args, self.guess, self.drop,
                              self.boundary, technique=self.technique)
      t2 = time.time()
      xout[3] = abs(xout[3]); xout[4] = abs(xout[4])
      print "After fitting:", time.ctime()
      print "Duration:", (t2-t1), "seconds"
      sys.stdout.flush()
      print "Best-fit parameters:"
      print xout
      self.xout = xout
      self.phistar_m = array(self.phistar_m)
      self.p_history = array(self.p_history)



if __name__ == "__main__":
   parfile = sys.argv[1]
   biv_obj = BivariateRL(parfile)
   biv_obj.mlfit()
   if len(sys.argv) > 2:  # save the state
      pname = sys.argv[2]
      biv_obj.pickle(pname)
   print "Done fitting."
