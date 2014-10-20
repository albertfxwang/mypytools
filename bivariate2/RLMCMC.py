#!/usr/bin/env python

import numpy as np
import emcee
import galaxy_samples as gs
import os, glob, sys, time
import cPickle

# MCMC stuff... mostly inherited from fit_lbg.py
today = time.strftime('%y%m%d')

class LBGSampleMCMC(gs.LBGSample):
   # A subclass of LBGSample that runs MCMC
   # For MCMC, one needs the best-fit parameters as the starting point 
   # (including phistar).
   def __init__(self, paramfile, filtername, nwalkers, **kwargs):
      super(LBGSampleMCMC, self).__init__(paramfile, filtername, **kwargs)
      self.sampler = None
      self.nwalkers = nwalkers

   def build_sampler(self, ndim=5):
      # use one thread for now... maybe figure out later how to use
      # multiprocessing
      # use emcee to perform MCMC sampling of the posterior probability distribution
      self.sampler = emcee.EnsembleSampler(self.nwalkers, ndim, self.mlfunc_mcmc, 
                                          threads=1)
      self.has_walkers = False

   def create_walkers(self, par0, par0_scale=[0.1,0.5,0.1,0.1,0.1]):
      alpha0 = np.random.uniform(par0[0]-par0_scale[0], 
                                 par0[0]+par0_scale[0], 
                                 size=self.nwalkers)
      mstar0 = np.random.uniform(par0[1]-par0_scale[1], 
                                 par0[1]+par0_scale[1], 
                                 size=self.nwalkers)
      logr00 = np.random.uniform(par0[2]-par0_scale[2], 
                              par0[2]+par0_scale[2], 
                              size=self.nwalkers)
      sigma0 = np.random.uniform(par0[3]-par0_scale[3], 
                                 par0[3]+par0_scale[3], 
                                 size=self.nwalkers)
      beta0 = np.random.uniform(par0[4]-par0_scale[4], 
                                par0[4]+par0_scale[4], 
                                size=self.nwalkers)
      par0 = np.concatenate([alpha0,mstar0,logr00,sigma0,beta0])
      par0 = par0.reshape(5, self.nwalkers).T
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
            progress_file = "mcmc_run_%s.txt" % today
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
            print "Iteration %d" % self.sampler.iterations
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
         # else:
         #    pos, prob, state, blobs = self.sampler.run_mcmc(self.mcmc_pos, 
                                                            # iterations)
         self.mcmc_pos = pos
         self.mcmc_prob = prob
         self.mcmc_state = state
