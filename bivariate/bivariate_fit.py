#!/usr/bin/env python
# bivariate_fit
#
# v0.5.0 H. Ferguson 1/5/05
# v1.0.0 K-H Huang   5/20/09
# ported from bivariate_multiple
# Last update: 07/06/09

__version__ = "3.0.0 22 May 2013"
__author__ = "Kuang-Han Huang, JHU/STScI"

import time
import sys
import numpy as N
import bivRL as bl
import bivRM as bm
try:
   from scipy import optimize
except ImportError:
   print "Cannot import scipy.optimize."
from transfunc import mlutil
from mlutil import grid_data
#import interlopers
#import zdist
#from pygoods import *
#import schechter
#import gauss
#import lognormal 
#import cPickle
#from pygoods.numprint import *

# Globals
transfer_file = ""
_first_call_dict = {}
kgrid_dict = {}
kernel_list_dict = {}
kernel_mask_dict = {}
kgrid = None
restlim = N.array([[-25.0, -15.0],
                 [ -2.0 ,   3.0]])
pixdx = N.array([ 0.02,  0.02])
limits = N.array([[23.0,27.0],[-0.5,2.0]])
datalimits=limits.copy()
goods_area = 160.  # in arcmin**2
udf_area = 11.     # in arcmin**2

##### Interloper fraction from photometric scatter, derived from 
##### simulations #####
interloper_frac_goods = {'v':[0.,0.,0.,0.,0.,0.0263,0.0614,0.0394,0.0468,
                              0.061,0.079],
                         'b':[0.,0.,0.,0.,0.0076,0.0032,0.0207,0.0216,
                              0.0276,0.0498,0.0728]}
interloper_mag_goods = N.arange(21., 26.5, 0.5)
interloper_frac_udf = {'v':[0.,0.,0.,0.,0.0988,0.0291,0.0546,0.0179,0.0206,
                            0.0394,0.0542],
                       'b':[0.,0.,0.0087,0.,0.0564,0.0242,0.0509,0.0503,
                            0.0537,0.0822,0.1164]}
interloper_mag_udf = N.arange(23.0, 28.5, 0.5)
interloper_frac_goods_all = {'v':0.0707, 'b':0.0550}
interloper_frac_udf_all = {'v':0.0377, 'b':0.076}
#############################################################################

### ===== For tests... ========================= ###

def rand_guess(plow=[-2.0,-30.,0.0,0.,0.1],phigh=[-0.5,-15.,2.0,2.0,1.0]):
   alpha = N.random.uniform(plow[0],phigh[0])
   Mstar = N.random.uniform(plow[1],phigh[1])
   logrpeak = N.random.uniform(plow[2],phigh[2])
   sigma = N.random.uniform(plow[3],phigh[3])
   beta = N.random.uniform(plow[4],phigh[4])
   return N.array([alpha,Mstar,logrpeak,sigma,beta])
    
def make_config(outname,guess,datafile='simdata.dat',technique="bfgs",
                image_scales=0.03,m0=25.0,transfer_fct='kernel5.p',fft=-1,
                min_comp=-1,limits=limits,dropdic='bdrops.p',
                datalimits=datalimits,pixdx=pixdx,verbose=1,
                kcorrfile='kcorr1700_i.p',flag=1,reerrcut=1.0,
                cullflags=[0,11,14]):
   f = open(outname,'w')
   # DATAFILE: GALFIT catalog of LBGs
   f.write('DATAFILE     ')
   f.write('%s ' % datafile)
   f.write('#\n')
   f.write('TECHNIQUE     %s\n' % technique)
   f.write('IMAGE_SCALES      %.3f\n' % image_scales)
   f.write('ALPHA      %f\n' % guess[0])
   f.write('MSTAR      %f\n' % guess[1])
   f.write('LOGRPEAK      %f\n' % guess[2])
   f.write('SIGMA      %f\n' % guess[3])
   f.write('M0         %f\n' % m0)
   f.write('BETA       %f\n' % guess[4])
   if len(transfer_fct)>0:
      f.write('TRANSFER_FCT      %s\n' % transfer_fct)
   if fft >= 0:
      f.write('FFT        %d\n' % fft)
   if min_comp >= 0:
      f.write('MIN_COMP          %f\n' % min_comp) 
   f.write('DROPDIC     %s\n' % dropdic)
   f.write('KCORRFILE   %s\n' % kcorrfile)
   f.write('MODMAG0     %f\n' % limits[0,0])
   f.write('MODMAG1     %f\n' % limits[0,1])
   f.write('MODLOGR0     %f\n' % limits[1,0])
   f.write('MODLOGR1     %f\n' % limits[1,1])
   f.write('DATAMAG0     %f\n' % datalimits[0,0])
   f.write('DATAMAG1     %f\n' % datalimits[0,1])
   f.write('DATALOGR0     %f\n' % datalimits[1,0])
   f.write('DATALOGR1     %f\n' % datalimits[1,1])
   f.write('PIXDX        %f   %f\n' % (pixdx[0],pixdx[1]))
   f.write('VERBOSE      %d\n' % verbose)
   f.write('FLAG         %d\n' % flag)
   f.write('REERRCUT     %d\n' % reerrcut)
   f.write('CULLFLAGS    ')
   for cf in cullflags:
      f.write('%d ' % cf)
   f.write('\n')
   f.flush()
   f.close()

### ===== The real fitting functions ================== ###                

class biv_fit(object):
   def loglikelihood(self,data,inputmodel,floor=1.e-50):
      """
      Return the -sum(log(likelhood))
      Arguments:
      data - data array of shape (ndata,ndim)
      inputmodel - 2-D model array (pdf normalized to same total number) 
      """
      if N.shape(data)[0]!=2:
         raise ValueError, "input data array has incorrect shape"
      #ndata = N.shape(data)[1] # number of data points
      model = N.maximum(inputmodel.model,floor)
      maglims = inputmodel.limits[0]
      logrlims = inputmodel.limits[1]
      pixdx = inputmodel.pixdx
      # Construct the pixel edges in both dimensions 
      mag_pix = N.arange(maglims[0], maglims[1]+pixdx[0], pixdx[0])
      logr_pix = N.arange(logrlims[0], logrlims[1]+pixdx[1], pixdx[1])
      datagrid = N.histogram2d(data[0], data[1], bins=[mag_pix, logr_pix])[0]
      #datagrid = mlutil.grid_data(data,inputmodel.limits,inputmodel.pixdx)

      if N.sum(datagrid.ravel()) == 0:
         #print "No data used; return 0."
         return 0
      model_l = N.compress(datagrid.ravel()>0, model.ravel())
      datagrid_l = N.compress(datagrid.ravel()>0, datagrid.ravel())
      logl = N.sum((datagrid_l*N.log(model_l)).ravel())
      #logl = N.sum((datagrid * N.log(model)).ravel())
      return -logl

   def mlfunc_RL_mcmc(self, parameters, *args):
      if self.sampler == None:
         raise ValueError, "Cannot call mlfunc_RL_mcmc if not running MCMC!!"
      p = self.mlfunc_RL(parameters, *args)
      return (-1.)*p[0], p[1]

   def calc_model_field(self, field, pars=None):
      ff = field
      if pars==None:
         pars = self.parameters
      if not self.add_interloper:
         for f in self.fields:
            self.mag_lolims[f] = 0.
      self.models[ff].bivariate_RL(pars, self.limits[ff], self.pixdx, 
                                  self.drop, ff, kgrid=self.tfkgrids[ff], 
                                  zdgrid=self.zdgrids[ff],
                                  mc=self.mc[ff], M0=self.M0,
                                  add_interloper=self.add_interloper,
                                  sd=self.sd[ff], intfrac=self.intfrac[ff], 
                                  mag_lolim=self.mag_lolims[ff],
                                  nproc_model=self.nproc_model,
                                  delta_z=self.delta_z)

   def mlfunc_RL(self, parameters, datasets, limits, pixdx, tfkgrids, zdgrids,
           M0, mc, z_mean, drop, add_interloper, sd, intfrac, mag_lolim,
           nproc_model, zdist_mag_udf, kgrid1d_udf, mag_array_udf, maglim_udf):
      """ 
      Compute the loglikelihood of a data given a model. Will be called 
      by the optimization routine to return loglikelihood.

      Arguments:
      parameters     -- parameters of the model
      datasets       -- dictionary of data; each entry contains a 2xN array
                        the keys are the name of the fields
      #restlimits    -- rest-frame limit of the model (obsolete)
      limits         -- dictionary of model limits for each field
      pixdx          -- pixel scale of the model (same for all fields)
      tfkgrids       -- dictionary of transfer function kernel grids
      #norm          -- value that the model normalizes to (default is 1.0)
      zdgrids        -- dictionary of P(z) kernel grids for each field
      M0             -- the reference abs mag (default to -21 mag); shoule be
                        a fixed value
      mc             -- an instance of mconvert for magnitude conversion 
                        from M to m
      z_mean         -- the mean redshift of the sample
      drop           -- the dropout band of the sample
      add_interloper -- whether to include the interloper component or not
                     (needs further revision to read in interloper fractions)
      sd             -- a sizedist class instance
      intfrac        -- a FITS table containing interloper fractions
      mag_lolim      -- the lower-limit in magnitude brighter than which 
                        do not add interloper model contributions
      zdist_mag_udf  -- the P(z) as a function of M1500 in UDF fainter than
                        limits['udf'][0][1] (to extend the magnitude limit
                        to where GALFIT is not robust anymore)
      """
      #verbose = 1
      fields = datasets.keys()
      # a list of all field names
      Nobj = {}
      for f in fields:
         Nobj[f] = N.shape(datasets[f])[1]
         # the total number of points in each field
      # a dictionary that holds the model for all fields
      for f in fields:
         #models[f] = bl.bivariate_RL_class()
         self.models[f].bivariate_RL(parameters, limits[f], pixdx, drop,
                                  f, kgrid=tfkgrids[f], zdgrid=zdgrids[f],
                                  mc=mc[f], M0=M0,
                                  add_interloper=add_interloper,
                                  sd=sd[f], intfrac=intfrac[f], 
                                  mag_lolim=mag_lolim[f],
                                  nproc_model=nproc_model,
                                  delta_z=self.delta_z)
         # Calculate the model for each field
   
      # calculate phistar considering all fields
      modelsum = 0.
      Ntot = 0.
      for f in fields:
         modelsum += N.sum(self.models[f].model.ravel()) * pixdx[0] * pixdx[1]
         Ntot += Nobj[f]
      phistar = float(Ntot) / float(modelsum)
      self.p_history += [parameters]
      self.phistar += [phistar]
      # Store the current phistar
      logl_tot = 0.
      # Total loglikelihood
      for f in fields:
         self.models[f].model = phistar * self.models[f].model
         logl_tot += self.loglikelihood(datasets[f], self.models[f], floor=0.)
      # logl is actualy (-1) * log(likelihood)!   

      if self.verbose:
         print "pars, log(L), phistar: ", parameters, logl_tot, phistar

      # Below is a kludge to extend i-dropout fitting down to ~30 mag or so
      # in UDF
      #if drop == 'f775w':
      #   modelext = bl.univariate_LF_class(limits['udf'], maglim_udf, 
      #                                     zdgrids['udf'].zarr)
      #   modelext.univariate_LF_all(parameters, mc['udf'], zdist_mag_udf)
      #   modelext.apply_transfer_1d(kgrid1d_udf)
      #   logl_ext = modelext.loglikelihood(mag_array_udf)
      #   logl_tot += logl_ext
      #   if verbose:
      #      print "logl_ext = %.2f" % logl_ext
      if self.sampler == None:
         return logl_tot
      else:
         return logl_tot, phistar

   def mlfunc_RM(self, parameters, datasets, limits, pixdx, tfkgrids,
           zdgrids, logM0, mc, z_mean, drop, add_interloper, sd, intfrac, 
           mag_lolim, M2M_kernel, nproc_model):
      """ 
      Compute the loglikelihood of a data given an RM distribution model. 
      Will be called by the optimization routine to return loglikelihood.

      Arguments:
      parameters     -- parameters of the model
      datasets       -- dictionary of data; each entry contains a 2xN array
                        the keys are the name of the fields
      limits         -- dictionary of model limits for each field
      pixdx          -- pixel scale of the model (same for all fields)
      tfkgrids       -- dictionary of transfer function kernel grids
      #norm          -- value that the model normalizes to (default is 1.0)
      zdgrids        -- dictionary of P(z) kernel grids for each field
      M0             -- the reference abs mag (default to -21 mag); shoule be
                        a fixed value
      mc             -- an instance of mconvert for magnitude conversion 
                        from M to m
      z_mean         -- the mean redshift of the sample
      drop           -- the dropout band of the sample
      add_interloper -- whether to include the interloper component or not
                     (needs further revision to read in interloper fractions)
      sd             -- a sizedist class instance
      intfrac        -- a FITS table containing interloper fractions
      mag_lolim      -- the lower-limit in magnitude brighter than which 
                        do not add interloper model contributions
      M2M_kernel     -- log-stellar mass to M1500 conversion kernel, based on
                        the relation M1500 = -logM_star + X, and X should
                        follow a distribution around X0.
      """
      verbose = 0
      #print "pars:", parameters
      fields = datasets.keys()
      # a list of all field names
      Nobj = {}
      for f in fields:
         Nobj[f] = N.shape(datasets[f])[1]
         # the total number of points in each field
      # a dictionary that holds the model for all fields
      for f in fields:
         #models[f] = bl.bivariate_RL_class()
         self.models[f].bivariate_RM(parameters, limits[f], pixdx,
                                     drop, f, kgrid=tfkgrids[f], 
                                     zdgrid=zdgrids[f], mc=mc[f], logM0=logM0,
                                     M2M_kernel=M2M_kernel, 
                                     add_interloper=add_interloper,
                                     sd=sd[f], intfrac=intfrac[f], 
                                     mag_lolim=mag_lolim[f],
                                     nproc_model=nproc_model)
         # Calculate the model for each field
   
      # calculate phistar_m considering all fields
      modelsum = 0.
      Ntot = 0.
      for f in fields:
         modelsum += N.sum(self.models[f].model.ravel()) * pixdx[0] * pixdx[1]
         Ntot += Nobj[f]
      phistar_m = float(Ntot) / float(modelsum)
      self.p_history += [parameters]
      self.phistar_m += [phistar_m]
      # Store the current phistar_m
      logl_tot = 0.
      # Total loglikelihood
      for f in fields:
         self.models[f].model = phistar_m * self.models[f].model
         logl_tot += self.loglikelihood(datasets[f], self.models[f], 
                                        floor=1.e-50)
      # logl is actualy (-1) * log(likelihood)!   

      if verbose:
         print "pars, log(L): ", parameters, logl_tot
      return logl_tot


   def bivariate_fit(self, mlfunc_args, guess, drop, boundary,
                        technique="l_bfgs_b", verbose=2, maxiter=200):
      """
      Cases for minimization technique:
      simplex:       not always robust
      anneal:        very slow
      bfgs:          good performance, but do not support constraints
      l_bfgs_b:      supports constraints
      test_pipe:     randomly generates an output (to test that everything 
                     before the minimization step is correct
      test_mlfunc:   returns the log-likelihood for the given values of the 
                     parameters
      Arguments:
      mlfunc_args:   arguments to mlfunc, in the exact order
      guess:         initial guess of the parameters
      drop:          name of the dropout band
      boundary:      parameter boundaries 
                     (an Nx2 array, N=number of parameters)
      technique:     minimization algorithm; possible options are simplex,
                     anneal, bfgs, l_bfgs_b
      verbose:       frequency to print intermediate steps of optimization
      maxiter:       max. number of iterations before reporting results
      """
      t1 = time.time()
   
      M0 = mlfunc_args[6]
   
      # ENTERING FITTING ALGORITHM

      # Simplex fit 
      if technique == "simplex":
         simplex_overrides = {
            "xtol":1.e-4,
            "ftol":1.e-6,
            "maxiter":1000,
            "full_output":1
            }

         sys.stdout.flush()
         print "ftol =", simplex_overrides["ftol"]
         (r,fopt,niters,funcalls,warnflag) = optimize.fmin(self.mlfunc,
            guess,args=mlfunc_args,**simplex_overrides)
         print "output array:", xopt
         print "alpha = %10.5f" % r[0]
         print "Mstar = %10.5f" % (r[1])
         print "rpeak = %10.5f pixels at M0 = " % (10.**r[2]), M0
         print "sigma = %10.5f" % (r[3])
         print "beta = %10.5f" % (r[4])
         print "-Log(L) = %10.2f" % (fopt)
         print "phistar = %.3e" % (self.phistar[-1])
         print "iter = %6d" % (niters)
         print r
         output = r
         t2 = time.time()
         print "Total fitting time: %.2f minutes" % ((t2-t1)/60.)
         return output

      #   Annealing fit
      #   Not updated... because it is too slow (130522)
      if technique == "anneal":
         if kgrid1 != None:
            print "uses kgrid1 in anneal fitting"
         else:
            print "no kgrid1"
         if kgrid2 != None:
            print "uses kgrid2 in anneal fitting"
         else:
            print "no kgrid2"

         anneal_args={}
         anneal_args['schedule']='fast'
         anneal_args['T0']=0.0001
         anneal_args['dwell']=1000
         anneal_args['boltzmann']=1.0
         anneal_args['full_output']=1
         anneal_args['maxeval']=10000
         anneal_args['maxiter']=1000
         anneal_args['feps']=1.e-30
         anneal_args['lower']=N.array([-2.5,21.,-1.,0.1,0.1]) 
         anneal_args['upper']=N.array([-0.5,28.,1.,1.0,1.0]) 
         anneal_args['args']=mlfunc_args
         (r,jmin,T,feval,iter,accept,retval) = optimize.anneal(self.mlfunc,
            guess,**anneal_args)
         print "output array:", xopt
         print "alpha = %10.5f" % (r[0])
         print "Mstar = %10.5f" % (r[1])
         print "rpeak = %10.5f pixels at M0 = " % (10.**r[2]), M0
         print "sigma = %10.5f" % (r[3]/N.log(10.))
         print "beta = %10.5f" % (r[4])
         print "-Log(L) = %10.2f" % (jmin)
         print "phistar = %.3e" % (self.phistar[-1])
         print "iter = %6d" % (iter)
         output = r
         t2 = time.time()
         print "Total fitting time: %.2f minutes" % ((t2-t1)/60.)
         return output

      # BFGS algorithm
      if technique == 'bfgs':

         print "maxiter =", maxiter
         gtol = 1.e-5
         print "gtol =",gtol
         sys.stdout.flush()
         all_output = optimize.fmin_bfgs( self.mlfunc, guess,
                      args=mlfunc_args, full_output=1, maxiter=maxiter,
                      gtol=gtol, retall=0)
         (xopt,fopt,gopt,Bopt,func_calls,grad_calls,warnflag) = all_output
         print "output array:", xopt
         print "alpha = %10.5f" % (xopt[0])
         print "Mstar = %10.5f" % (xopt[1])
         print "rpeak = %10.5f pixels at M0 = " % (10.**xopt[2]), M0
         print "sigma = %10.5f" % (xopt[3])
         print "beta = %10.5f" % (xopt[4])
         print "-Log(L) = %10.2f" % (fopt)
         print "phistar = %.3e" % (self.phistar[-1])
         #print "func_calls = %6d" % (func_calls)
         #print "grad_calls = %6d" % (grad_calls)
         output = xopt
         t2 = time.time()
         print "Total fitting time: %.2f minutes" % ((t2-t1)/60.)
         return output

      if technique == 'l_bfgs_b':
         # parameter boundaries
         print "boundary:", boundary

         factr = 1.e7
         # 1e12 for low accuracy; 1.e7 for moderate accuracy; 
         # 10.0 for high accuracy
         print "factor: %.2e" % factr
         sys.stdout.flush()

         # Do optimization
         (xopt,fopt,d_fit) = optimize.fmin_l_bfgs_b(self.mlfunc,
            guess,args=mlfunc_args,fprime=None,approx_grad=True,
            factr=factr,epsilon=1.e-5,bounds=boundary,iprint=verbose)

         # Print results
         print "output array:", xopt
         print "alpha = %.5f" % xopt[0]
         print "Mstar = %10.5f" % xopt[1]
         print "rpeak = %10.5f pixels at M0 = " % (10.**xopt[2]), M0
         print "sigma = %10.5f" % (xopt[3])
         print "beta = %10.5f" % (xopt[4])
         print "-Log(L) = %10.2f" % (fopt)
         print "phistar = %.3e" % (self.phistar[-1])
         print "fitting information:", d_fit
         output = xopt
         t2 = time.time()
         print "Total fitting time: %.2f minutes" % ((t2-t1)/60.)
         return output

      # Not updated (130522)
      if technique == "test_pipe":
         print "boundary", boundary
         print "guess", guess
         xout = N.zeros(len(guess))
         for i in range(len(guess)):
            if None in boundary[i]:
               xout[i] = N.random.normal(guess[i],1.0)
            else:
               xout[i] = N.random.uniform(boundary[i][0],boundary[i][1])
         print "xout:",xout
         return xout

      # Not updated (130522)
      if technique == "test_mlfunc":
         value = self.mlfunc(guess, *mlfunc_args)
         print value
