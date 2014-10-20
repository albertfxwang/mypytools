#!/usr/bin/env python

# Under bivariate2
# For Maximum-likelihood fitting...
__version__ = "1.0"

import numpy as np
from scipy import optimize
import galaxy_samples
import sys, time

# available minimization algorithms
available_techniques = ['simplex', 'anneal', 'bfgs', 'l_bfgs_b']

class LBGSampleFitting(galaxy_samples.LBGSample):
   # a subclass for Maximum-likelihood fitting
   def print_output(self):
      # print "output array:", xopt
      assert len(self.output)==6, "Something wrong with self.output... it has only %d elements." % len(self.output)
      r = self.output
      print "alpha = %10.4f" % r[0]
      print "Mstar = %10.5f" % (r[1])
      print "logR0 = %10.5f (R0 in arcsec)" % (r[2])
      print "sigma = %10.5f" % (r[3])
      print "beta = %10.5f" % (r[4])
      print "phistar = %10.3e" % (r[5])
      print "-Log(L) = %10.2f" % (self.maxlikelihood)
      m, s = divmod(self.dt, 60.)
      h, m = divmod(m, 60.)
      print "Total fitting time: %d hour %d minutes %.2f seconds" % (h, m, s)

   def simplex_fitting(self):
      # Simplex minimization; reasonably fast
      print "Minimization with Simplex..."
      simplex_overrides = {
            "xtol":1.e-4,
            "ftol":1.e-6,
            "maxiter":1000,
            "full_output":1
            }
      sys.stdout.flush()
      print "ftol =", simplex_overrides["ftol"]
      (r,fopt,niters,funcalls,warnflag) = optimize.fmin(self.mlfunc, self.guess,
         **simplex_overrides)
      self.output = np.zeros(6)
      self.output[:5] = r
      self.output[5] = self.phistar
      self.maxlikelihood = fopt  # it's actually minimum -logl

   def anneal_fitting(self):
      # Annealing fit
      # Not updated... because it is too slow (130522)
      print "Minimization with Annealing..."
      anneal_args={}
      anneal_args['schedule'] = 'fast'
      anneal_args['T0'] = 0.0001
      anneal_args['dwell'] = 1000
      anneal_args['boltzmann'] = 1.0
      anneal_args['full_output'] = 1
      anneal_args['maxeval'] = 10000
      anneal_args['maxiter'] = 1000
      anneal_args['feps'] = 1.e-30
      anneal_args['lower'] = N.array([-2.5,21.,-1.,0.1,0.1]) 
      anneal_args['upper'] = N.array([-0.5,28.,1.,1.0,1.0]) 
      anneal_args['args'] = mlfunc_args
      (r, jmin, T, feval, niter, accept, retval) = optimize.anneal(self.mlfunc,
         self.guess, **anneal_args)
      self.output = np.zeros(6)
      self.output[:5] = r
      self.output[5] = self.phistar
      self.maxlikelihood = jmin

   def bfgs_fitting(self, maxiter=200):
      # BFGS algorithm
      print "Minimization with BFGS..."
      print "maxiter =", maxiter
      gtol = 1.e-5
      print "gtol =",gtol
      sys.stdout.flush()
      all_output = optimize.fmin_bfgs( self.mlfunc, self.guess, full_output=1, 
                                      maxiter=maxiter, gtol=gtol, retall=0)
      (xopt,fopt,gopt,Bopt,func_calls,grad_calls,warnflag) = all_output
      self.output = np.zeros(6)
      self.output[:5] = xopt
      self.output[5] = self.phistar
      self.maxlikelihood = fopt
      
   def l_bfgs_b_fitting(self, factr=1.e7, boundary=None, verbose=2):
      # L_BFGS_B algorithm... able to take boundaries, too, but not implemented yet
      # if boundary != None:
      #    raise NotImplementedError
      # factr controls the termination criteria.
      # 1e12 for low accuracy; 1.e7 for moderate accuracy; 
      # 10.0 for high accuracy
      print "factor: %.2e" % factr
      sys.stdout.flush()
      # Do optimization
      ### Add boundary option later...
      (xopt,fopt,d_fit) = optimize.fmin_l_bfgs_b(self.mlfunc, self.guess,
                           fprime=None, approx_grad=True, factr=factr,
                           epsilon=1.e-5, iprint=verbose, bounds=boundary)
      self.output = np.zeros(6)
      self.output[:5] = xopt
      self.output[5] = self.phistar
      self.maxlikelihood = fopt
   
   def max_likelihood_fitting(self, technique=None, **kwargs):
      """
      The driver method for max-likelihood fitting.
      Keyword arguments are fed to individual fitting methods.
      """
      # Print all fitting parameters
      for key in self.p.keys():
         print key, self.p[key]
      if technique == None:
         technique = self.technique
      if technique not in available_techniques:
         raise NotImplementedError, "Technique %s not yet implemented." % technique
      t1 = time.time()
      if technique == 'l_bfgs_b':
         self.l_bfgs_b_fitting(**kwargs)
      elif technique == 'bfgs':
         self.bfgs_fitting(**kwargs)
      elif technique == 'simplex':
         self.simplex_fitting(self)
      elif technique == 'anneal':
         self.anneal_fitting(self)
      t2 = time.time()
      self.dt = t2 - t1
      self.print_output()

if __name__ == "__main__":
   paramfile = sys.argv[1]
   filtername = sys.argv[2]
   ## does not allow extra keyword arguments at the moment...
   LBG = LBGSampleFitting(paramfile, filtername)
   LBG.max_likelihood_fitting()

      