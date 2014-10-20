#!/usr/bin/env python
# version 2 -- use emcee (MCMC Hammer) instead of pymc; much more intuitive (and faster??)

from numpy import *
#import pymc
import emcee
import os, glob, sys, time
import bivariate_fit as bf
#import mlutil
from pygoods import parseconfig
import bivRL as bl
import cPickle

"""
Use emcee to run Markov Chain Monte Carlo in order to sample the posterior 
probability distribution of the parameters in the bivariate model. 5 parameters
are alpha, mstar, r0, sigma, and beta. Start with the best-fit parameter
values found by bivariate_fit.py (using one of the minimization algorithms
provided by scipy, default is L_BFGS_B). Set proper burn-in period and obtain
posterior distributions of each parameter.
"""

def build_sampler(mlfunc_args, nwalkers, ndim=5, nthreads=1):
	# use emcee to perform MCMC sampling of the posterior probability distribution
	sampler = emcee.EnsembleSampler(nwalkers, ndim, bf.mlfunc_RL, 
                                   args=mlfunc_args, threads=nthreads)
	return sampler
	

def readpar(parfile):
   c = parseconfig(parfile)
   for k in c.keys():
      print k, c[k]
   limits1 = array(c['LIMITS1'])
   c['LIMITS1'] = limits1.reshape((2,2))
   limits2 = array(c['LIMITS2'])
   c['LIMITS2'] = limits2.reshape((2,2))
   if c['MCFILE'] == "None": mcfile = None
   if c['ZDGRIDFILE1'] == "None": c['ZDGRIDFILE1'] = None
   if c['ZDGRIDFILE2'] == "Nene": c['ZDGRIDFILE2'] = None
   bound_alpha = c['BOUND_ALPHA']
   c['BOUND_ALPHA'] = replaceNone(bound_alpha)
   bound_Mstar = c['BOUND_MSTAR']
   c['BOUND_STAR'] = replaceNone(bound_Mstar)
   bound_r0 = c['BOUND_R0']
   c['BOUND_R0'] = replaceNone(bound_r0)
   bound_sigma = c['BOUND_SIGMA']
   c['BOUND_SIGMA'] = replaceNone(bound_sigma)
   bound_beta = c['BOUND_BETA']
   c['BOUND_BETA'] = replaceNone(bound_beta)
   figname = c['FIGNAME']
   plot = c['PLOT']
   if c['PLOT'] == "yes": c['PLOT'] = True
   else: c['PLOT'] = False
   #if (c['FIT2'] == "True") | (c['FIT2'] == "yes"):
   #   c['FIT2'] = True
   #else:
   #   c['FIT2'] = False 
   c['CULLFLAGS'] = array(c['CULLFLAGS'])
   return c

def bivmodel_MCMC_pymc(par, parfile, picklename, db='pickle', zlo=-1):
   """
   par: best-fit parameter found by bivariate_fit.py
   parfile: parameter file name of fitting run
   """
   c = readpar(parfile)
   
   # define 5 independent variables
   alpha = pymc.Uniform(name='alpha', lower=-4.0, upper=0.0)
   alpha.set_value(par[0])
   Mstar = pymc.Uniform(name='Mstar', lower=-25.0, upper=-15.0)
   Mstar.set_value(par[1])
   logr0 = pymc.Uniform(name='logr0', lower=0.0, upper=5.0)
   logr0.set_value(par[2])
   sigma = pymc.Uniform(name='sigma', lower=0.0, upper=5.0)
   sigma.set_value(par[3])
   beta = pymc.Uniform(name='beta', lower=0.0, upper=10.0)
   beta.set_value(par[4])
   # define data
   #if c['DROP'] == 'b': zlo = 3.0
   #else: zlo = 4.0
   mag1,re1,crit1 = cleandata(c['LBGCAT1'],
                              reerr_ratiolim=c['REERR_RATIOLIM1'],
                              chisqnulim=c['CHISQNULIM1'],
                              cullflags=c['CULLFLAGS'],
                              magautolim=c['MAGAUTOLIM1'],
                              limits=c['LIMITS1'],
                              drop=c['DROP'], zlo=zlo)
   data1 = array([mag1,log10(re1)]); print "shape(data1)", shape(data1)
   n1 = shape(data1)[1] # n1 is the number of dropouts in the ACS dataset
   mag2,re2,crit2 = cleandata(c['LBGCAT2'],
                              reerr_ratiolim=c['REERR_RATIOLIM2'],
                              chisqnulim=c['CHISQNULIM2'],
                              cullflags=c['CULLFLAGS'],
                              magautolim=c['MAGAUTOLIM2'],
                              limits=c['LIMITS2'],
                              drop=c['DROP'], zlo=zlo)
   data2 = array([mag2,log10(re2)]); print "shape(data2)", shape(data2)
   
   kgrid1 = None; kgrid2 = None
   if c['KERNEL1'] != "None":
      f = open(c['KERNEL1'], 'rb')
      kgrid1 = cPickle.load(f)
      f.close()
   if c['KERNEL2'] != "None":
      f = open(c['KERNEL2'], 'rb')
      kgrid2 = cPickle.load(f)
      f.close()


   pixdx = array([0.02, 0.02])
   mc = bl.mconvert(c['MCFILE'])
   if c['DROP'] == 'b':
      meankcorr = mc(4.0)
   else:
      meankcorr = mc(5.0)
   mlfunc_arg = [c['LIMITS1'], c['LIMITS2'], pixdx, kgrid1, kgrid2,\
                 c['WEIGHT1'], c['WEIGHT2'], c['M0'], c['ZDGRIDFILE1'],\
                 c['ZDGRIDFILE2'], c['MCFILE'], 1, meankcorr, 0, -1., 'phistar',\
                 c['DROP']]

   # define data
   logl0 = bf.mlfunc_RL(par, data1, data2, *mlfunc_arg)
   value = concatenate([data1, data2], axis=1)
   print shape(value) #, shape(value[:,0:n1]), shape(value[:,n1:])
   @pymc.stochastic(observed=True)
   def D(alpha=alpha, Mstar=Mstar, logr0=logr0, sigma=sigma, beta=beta,
      value=value):
      def logp(value, alpha, Mstar, logr0, sigma, beta):
         logl = bf.mlfunc(array([alpha, Mstar, logr0, sigma, beta]),
            value[:,:n1], value[:,n1:], *mlfunc_arg)
         return (-1.) * (logl-logl0)
         #return (-1.) * logl

   if db == 'pickle':
      kwarg = {'db':db, 'dbname':picklename}
   else:
      kwarg = {'db':db}
   M = pymc.MCMC({'alpha':alpha, 'Mstar':Mstar, 'logr0':logr0, 'sigma':sigma, 'beta':beta,\
                  'D':D}, **kwarg)
   M.use_step_method(pymc.Metropolis, M.Mstar, scale=1., proposal_sd=0.1,
      proposal_distribution="Normal")
   M.use_step_method(pymc.Metropolis, M.alpha, scale=1., proposal_sd=0.1,
      proposal_distribution="Normal")
   M.use_step_method(pymc.Metropolis, M.logr0, scale=1., proposal_sd=0.1,
      proposal_distribution="Normal")
   M.use_step_method(pymc.Metropolis, M.sigma, scale=1., proposal_sd=0.1,
      proposal_distribution="Normal")
   M.use_step_method(pymc.Metropolis, M.beta,  scale=1., proposal_sd=0.1,
      proposal_distribution="Normal")

   return M, alpha, Mstar, logr0, sigma, beta, D


def sample_pymc(M, iter=100, burn=0, thin=1):
   t1 = time.time()
   M.sample(iter=iter, burn=burn, thin=thin, tune_throughout=True, tune_interval=20)
   t2 = time.time()
   dt = (t2-t1)/60.
   print "%d iterations took %.1f minutes" % (iter, dt)
   M.db.close()  # write db to designated backend
   return M


if __name__ == "__main__":
   c = parseconfig(sys.argv[1])
   hostname = os.getenv('HOSTNAME')
   par = array([c['ALPHA'], c['MSTAR'], c['LOGR0'], c['SIGMA'], c['BETA']])
   print "starting parameters:",par
   print "iterations:", c['ITER']
   if len(glob.glob(c['PICKLENAME'])):
      raise ValueError, "%s already exists" % c['PICKLENAME']
   M, alpha, Mstar, logr0, sigma, beta, D = bivmodel_MCMC(par, c['PARFILE'], c['PICKLENAME'])
   M = sample(M, iter=c['ITER'], burn=c['BURN'], thin=c['THIN'])
   #gmail.gmail('astrokuang@gmail.com','kuanghan@pha.jhu.edu',\
   #   '%s finished on %s' % (sys.argv[1], hostname))

