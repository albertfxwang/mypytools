#!/usr/bin/env python

from numpy import *
import scipy
from scipy import optimize, integrate, interpolate

## Author: Kuang-Han Huang
## Last updated: 2014/08/29

# use the Bayesian formulation of Press 1997 to calculate the posterior distribution
# of P(z0|D), where z0 is the quantity we want to measure, and D is the set of data and errors
# reported by different measurements
# We use it to generate posterior distribution of photo-z in our first application

class pdf:
   # a base class for probability density function
   def __init__(self, p, x):
      assert len(p) == len(x), 'Length of p (%d) and x (%d) are not compatible.' % (len(p), len(x))
      self.pdf = p
      self.x = x
      self.f = interpolate.interp1d(self.x, self.pdf)
      
   def __getitem__(self, xnew):
      # returns the probability density at xnew
      # use scipy to interpolate
      xnew = maximum(xnew, min(self.x))
      xnew = minimum(xnew, max(self.x))
      return self.f(xnew)
      
   def __call__(self, y):
      return self.__getitem__(y)
      
   def __len__(self):
      return len(self.pdf)

def gauss(x, m, sigma):
   sigma2 = sigma**2
   y = exp(-(x - m)**2/(2. * sigma2)) / sqrt(2. * pi * sigma2)
   y = (1. / sqrt(2 * pi * sigma2)) * y
   return y

def conf_inter_sym(p, val, int_guess=0.1, frac=0.68):
   # Given a PDF p, calculate a confidence interval around the desired value 
   # val. The total probability contained within the interval is frac.
   # p -- the input PDF function (need not be normalized); should be a pdf class instance
   # val -- the central value around which we want the confidence interval
   # int_guess -- the initial guess of the width of the confidence interval
   #    so that val +- int_guess brackets frac*(total probability)
   # frac -- the fraction of total probability the confidence interval should bracket
   #p.pdf = p.pdf / sum(p.pdf)  # normalize p so that sum(p.pdf) = 1
   ptot = integrate.quad(p, p.x[0], p.x[-1])[0]
   def mfunc(conf_inter):
      xlo = val - conf_inter
      xhi = val + conf_inter
      result = integrate.romberg(p, xlo, xhi)
      #print xlo, xhi, result
      psum = result
      #print psum
      #if result[1]/result[0] >= 1.e-3:
      #   print "Warning: fractional error of integration is too large! %.2f" % (result[1]/result[0])
      return abs(psum - frac*ptot)  
      # the difference between the bracketed probability and the desired fraction
   best_int = optimize.fmin(mfunc, int_guess, disp=0, full_output=0)[0]
   return best_int
   
   
class uniform_prior(pdf):
   # a class for uniform (uninformative) prior (not normalized)
   def __init__(self, x):
      pdf.__init__(self, ones(len(x)), x)
      
# class bayesian_data:
class BayesianPz(object):
   # a class that calculates the posterior distribution
   # assuming that each individual photo-z code returns a best photo-z with
   # a Gaussian error sigma...
   def __init__(self, data, zarray, sigma=None, pfloor=0.7, alpha=2.1):
      self.data = array(data) 
      # self.data could be an array of N photo-z estimates, or an N by Nz
      # array with each row being the individual P(z) estimate; Nz is the 
      # number of redshift steps in P(z).
      # If data is a 1-D, length N array, sigma needs to be a length N array
      # specifying the Gaussian errors around the best photo-z estimates.
      # First, if the input is N photo-z estimates and photo-z error estimates
      if self.data.ndim == 1:
         # assume that errors are Gaussian
         self.sigma = array(sigma)
         self.datamode = 'gaussian'
      else:
         # if the input is a collection of P(z)
         self.datamode = 'Pz'
         # re-format self.data into a list of pdf instances
         dataPz = []
         for i in range(len(data)):
            dataPz += [pdf(data[i], zarray)]
         self.data = dataPz
      self.N = len(data)  # number of independent measurements
      self.zarray = zarray  # define the outpu redshift steps
      self.dz = (self.zarray[1] - self.zarray[0])  
      self.zmax = self.zarray.max()
      self.alpha = alpha
      # assume the redshift sampling is uniform...
      # Set a flat redshift prior
      self.zPrior = pdf(ones(len(self.zarray)), self.zarray)      
      # Prior on z; usually flat between the allowed redshift range 
      # Prior on p or f_good or (1 - f_bad)
      # self.S = self.calc_S()
      # self.S is either the normalization if P(z)_bad is flat, or the 
      # inflated Gaussian error if P(z)_bad is gaussian.
      self.pRange = linspace(0, 1.0, 101)
      # set a prior of p (= f_good = 1 - f_bad); prior is 1 for p >= pfloor;
      # otherwise prior of p = 0
      gPrior = zeros(len(self.pRange))
      gPrior[self.pRange >= pfloor] = 1.0
      self.goodPrior = pdf(gPrior, self.pRange)  
      # the sampling of p or f_good
      self.dp = (self.pRange[1] - self.pRange[0])
      # Now initialize self.Pgood and self.Pbad
      # Both are a M by N matrix, where M is the number of photo-z estimates
      # and N is the number of redshift steps
      self.Pgood = zeros((len(self.data), len(self.zarray)))
      # Use flat P(z) for bad photo-z estimates... could try other forms 
      # in the future
      self.Pbad = ones(self.Pgood.shape) / self.zmax 
      # set up self.Pgood for each photo-z estimate 
      for i in range(len(self.data)):
         if self.datamode == 'gaussian':
            self.Pgood[i] = gauss(self.zarray, self.data[i], self.sigma[i])
         else:
            self.Pgood[i] = self.data[i](self.zarray)

   # def calc_S_gauss(self, factor=100.):
   # def calc_S(self, factor=100.):
   #    if self.datamode == 'gaussian':
   #       S = max(self.sigma) * factor  # the assumed error for bad measurements
   #    else:
   #       S = 1.0 / max(self.zarray)
   #    return S

   # def setGoodPrior(self, zarrPrior):
   #    # self.prior_p = pdf(zarrPrior, self.p_range)
   #    self.goodPrior = pdf(zarrPrior, self.p_range)
   #    # set prior on the probability that a given measurement is correct
   #    # In other words, p = f_good = 1 - f_bad
         
   # def gauss_prob(self, zi, z0, sigma):
   #    prob = exp(-1 * (zi - z0)**2 / (2.0 * sigma**2))
   #    # normalize the probability 
   #    prob = (1. / sqrt(2 * pi * sigma**2))
   #    return 
         
   # def calc_PBG_gauss(self, p, z0, pbad='flat'):
   def calc_Pgood_Pbad(self, p):
      # calculate the term prod(p*PG+(1-p)*PB) at all redshifts
      # pbad specifies what to use as PB ('flat' for a flat P(z) or 
      # 'gaussian' for a Gaussian with an inflated error).
      # if self.datamode == 'gaussian':
      #    Pgood_i = self.gauss_prob(self.data, zi, self.sigma)  # an array of length self.N
      #    Pbad_i = self.gauss_prob(self.data, zi, self.S) 
      # else:
      #    # take the input P(z) values at z0
      #    Pgood_i = array([P(zi) for P in self.data])
      #    Pbad_i = ones(self.N, 'float') * self.S / self.zmax
      #    # a flat PB at z = zi
      X = p * self.Pgood + (1. - p) * self.Pbad
      # returns a length N array, where N = len(self.zarray)
      return power(prod(X, 0), 1. / self.alpha)  
      
   # def calc_posterior_z0(self):
   #    # calculate the posterior probability density P(z=zi|D)
   #    # P(z=zi|D) = P(zi)*sum_{p,v}(P(p)*prod(p*PG+(1-p)*PB))
   #    # the term within prod() is calculated by self.calc_Pgood_Pbad
   #    # marginalize over all p and possible v
   #    # P1 = array([self.goodPrior(p) * self.calc_Pgood_Pbad(p, alpha=alpha) for p in self.pRange])
   #    P1 = reduce(lambda p: self.calc_Pgood_Pbad(p),
   #                self.pRange)
   #    return sum(P1), P1
   
   # def calc_posterior_p0(self, p0):
   #    # calculate the posterior distribution P(p0|D)
   #    # using self.zPrior as the prior on z
   #    # for this, use a flat prior on p: Prior(p) = 1 for all p; could experiment with other priors?
   #    # P2 = array([self.zPrior(z) * self.calc_Pgood_Pbad(p0, z) for z in self.zarray])
   #    P2 = reduce(lambda z: self.zPrior(z) * self.calc_Pgood_Pbad(p0, z), 
   #                self.zarray)
   #    return sum(P2), P2
      
      
   def calcPosteriorPz(self, alpha=2.1):
      # calculate the posterior distribution P(z|D); THIS IS WHAT WE WANT!!

      # Pz_post = array([self.calc_posterior_z0(z)[0] for z in self.zarray])
      Pz_map = map(lambda p: self.goodPrior(p) * self.calc_Pgood_Pbad(p), 
                   self.pRange)
      Pz_post = reduce(lambda p1,p2: p1 + p2, Pz_map)
      # Pz_post = zeros(len(self.zarray))
      # for p in self.pRange:
      #    Pz_post = Pz_post + self.calc_Pgood_Pbad(p)
      Pz_post = array(Pz_post)
      # normalize Pz_post to integrate to 1
      Pz_post = Pz_post  / (Pz_post.sum() * self.dz)
      self.posteriorPz = pdf(Pz_post, self.zarray)
      zPeak = self.zarray[argmax(Pz_post)]
      return self.posteriorPz, zPeak
   
   # def calcPosteriorPgood(self):
   #    # calculate P(p|D)
   #    # using self.prior_z as prior on z, which is basically the sum of P(z) from each measurement
   #    # Pgood_post = array([self.calc_posterior_p0(p)[0] for p in self.pRange])
   #    # normalize Pgood_post
   #    Pgood_post = Pgood_post / (Pgood_post.sum() * self.dp)
   #    self.posteriorPgood = Pgood_post
   #    pgoodPeak = self.pRange[argmax(Pgood_post)]  # the most likely value of p_good
   #    return self.posteriorPgood, pgoodPeak
      