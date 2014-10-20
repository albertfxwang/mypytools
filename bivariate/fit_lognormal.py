#!/usr/bin/env python

from numpy import *
import scipy
from scipy import optimize
import fit_lbg as fl
import bivRL as bl


# fit a lognormal distribution ** with natural log ** to the raw distribution


def logl_lognormal((r0, lnsigma), re):
   # returns the log-likelihood
   f = 1./(re * sqrt(2.*pi) * lnsigma) * exp(-1.*(log(re/r0)**2)/(2.*lnsigma**2))
   logl = log(f)
   return -1.*sum(logl)


#def fit_lognormal(drop='b', r0_guess=5.0, lnsigma_guess=0.5):
def fit_lognormal(re_array, r0_guess=5.0, lnsigma_guess=0.5):
   """
   Fit a log-normal function to the distribution of Re.
   """
   # r0_guess is in pixels
   guess = array([r0_guess, lnsigma_guess])
   print guess
   xout = optimize.fmin(logl_lognormal, guess, args=[re_array])
   return xout
   
def fit_lognormal_se(c1, c2, r0_guess=5.0, lnsigma_guess=0.5):
   re = concatenate((c1.flux_radius_1, c2.flux_radius_1))
   guess = array([r0_guess, lnsigma_guess])
   xout = optimize.fmin(logl_lognormal, guess, args=[re])
   return xout
      

def lognormal(r0, lnsigma, rarr):
   f = 1./(rarr * sqrt(2.*pi) * lnsigma) * exp(-1.*(log(rarr/r0)**2)/(2.*lnsigma**2))
   return f


