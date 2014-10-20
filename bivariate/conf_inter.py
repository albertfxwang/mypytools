#!/usr/bin/env python

from numpy import *
import scipy
from scipy import optimize
#from plot_mcmc import parb, parv
from mypytools import cosmoclass

cc = cosmoclass.cosmoclass(H0=70.2, omega_m=0.275, omega_l=0.725)

# from an array of parameter values and a best-fit value, find the shortest 68% (or any other
# percentile) interval around the best-fit value


def cdf(sortedN, x):
   # return the cumulative distribution function at x
   # does the array have to be sorted? maybe not...
   n_tot = len(sortedN)
   n_below = sum(sortedN <= x)
   return float(n_below)/float(n_tot)

def brac_percent(x1, x2, sortedN):
   y = cdf(sortedN, x2) - cdf(sortedN, x1)
   #print y
   return y

def dp(xhi, xlo, sortedN, p):
   return abs(p - brac_percent(xlo, xhi, sortedN))

def bracket(sortedN, xlo, p):
   # for a given value xlo and a sorted array of value sortedN, find the end point
   # xhi such that [xlo, xhi] bracket p percent of total number of points
   
   #dx = (max(sortedN)-min(sortedN)) / 2.
   dx = (xlo - min(sortedN)) / 5. # a rough step size
   #epsilon = (max(sortedN)-min(sortedN)) / len(sortedN)
   fout = optimize.fmin(dp, xlo+dx, args=(xlo, sortedN, p), disp=0,\
      full_output=1, xtol=1.e-4, ftol=1.e-4)
   #fout = optimize.brute(dp, ((xlo, max(sortedN)),), full_output=1)
   #print fout[:2]
   xopt = fout[0]
   return xopt


def conf_inter(sortedN, xcent, x0, p):
   # find the shortest interval around xcent that brackets a fraction p of the total 
   # probability
   # x0 is the initial guess of the lower endpoint of this interval
   # returns x0 and x1 so that the shortest interval is [x0, x1], which should bracket xcent
   if x0 > xcent: raise ValueError, "x0 cannot be larger than xcent"
   if x0 < min(sortedN): raise ValueError, "xlo can't be smaller than min(sortedN)"
   def len_inter(xlo, sortedN1, p1):
      xhi = bracket(sortedN1, xlo, p1)
      # find the higher endpoint that bracket probability p1
      if (xhi-xcent)*(xlo-xcent)>0:
         # if both xlo and xhi are on the same side of xcent, reject the result by returning the maximum distance
         return max(sortedN)-min(sortedN)
      else:
         return abs(xhi - xlo)  # just return the interval between xhi and xlo1
      #return (xhi-xcent)**2+(xcent-xlo)**2 # return the squared distances between xcent and both end points
      #return abs(xhi+xlo-2*xcent)  # return the distance difference between xcent to both end points

   fout = optimize.fmin(len_inter, x0, args=(sortedN, p), disp=0, full_output=1)
   #print fout
   xopt = fout[0][0]  # the lower endpoint
   #fopt = fout[1]  # what the value of len_inter is at xopt
   xopt_hi = bracket(sortedN, xopt, p)
   #print dp
   #xminus = xopt - xcent
   #xplus = xopt + fopt - xcent
   return array([xopt, xopt_hi])
   
      
def conf_inter_bdrops(c, bestpars, p=0.68):
   # alpha
   sorted_alpha = sort(c.alpha)
   x = conf_inter(sorted_alpha, bestpars[0], -1.8, p)
   print "alpha:", bestpars[0], x-bestpars[0], '\n'
   # mstar
   sorted_mstar = sort(c.mstar)
   x = conf_inter(sorted_mstar, bestpars[1], -20.8, p)
   print "mstar:", bestpars[1], x-bestpars[1], '\n'
   # phistar
   sorted_phistar = sort(c.phistar)
   x = conf_inter(sorted_phistar, bestpars[2], 1.2e-3, p)
   print "1000.*phistar:", bestpars[2]*1000., (x-bestpars[2])*1000., '\n'
   # logr0
   sorted_phistar = sort(c.logr0)
   x = conf_inter(sorted_phistar, bestpars[3], 0.8, p)
   print "logr0:", bestpars[3], (x-bestpars[3]), '\n'
   # logr0 in arcsec
   x_as = 10.**array(x) * 0.03
   r0_as = 10.**bestpars[3] * 0.03
   print "R0 [arcsec]:", r0_as, (x_as-r0_as), '\n'
   # logr0 in kpc
   x0_kpc = cc.adsize(4.0, x_as[0], unit='kpc')
   x1_kpc = cc.adsize(4.0, x_as[1], unit='kpc')
   r0_kpc = cc.adsize(4.0, r0_as, unit='kpc')
   print "R0 [kpc]:", r0_kpc, x0_kpc-r0_kpc, x1_kpc-r0_kpc, '\n'
   # sigma
   sorted_sigma = sort(c.sigma)
   x = conf_inter(sorted_sigma, bestpars[4], 0.7, p)
   print "sigma:", bestpars[4], (x-bestpars[4]), '\n'
   # beta
   sorted_beta = sort(c.beta)
   x = conf_inter(sorted_beta, bestpars[5], 0.2, p)
   print "beta:", bestpars[5], (x-bestpars[5]), '\n'
   print x


def conf_inter_vdrops(c, bestpars, p=0.68):
   # alpha
   sorted_alpha = sort(c.alpha)
   x = conf_inter(sorted_alpha, bestpars[0], -1.8, p)
   print "alpha:", bestpars[0], x-bestpars[0], '\n'
   # mstar
   sorted_mstar = sort(c.mstar)
   x = conf_inter(sorted_mstar, bestpars[1], -20.6, p)
   print "mstar:", bestpars[1], x-bestpars[1], '\n'
   # phistar
   sorted_phistar = sort(c.phistar)
   x = conf_inter(sorted_phistar, bestpars[2], 1.0e-3, p)
   print "1000.*phistar:", bestpars[2]*1000., (x-bestpars[2])*1000., '\n'
   # logr0
   sorted_phistar = sort(c.logr0)
   x = conf_inter(sorted_phistar, bestpars[3], 0.8, p)
   print "logr0:", bestpars[3], (x-bestpars[3]), '\n'
   # logr0 in arcsec
   x_as = 10.**array(x) * 0.03
   r0_as = 10.**bestpars[3] * 0.03
   print "R0 [arcsec]:", r0_as, (x_as-r0_as), '\n'
   # logr0 in kpc
   x0_kpc = cc.adsize(5.0, x_as[0], unit='kpc')
   x1_kpc = cc.adsize(5.0, x_as[1], unit='kpc')
   r0_kpc = cc.adsize(5.0, r0_as, unit='kpc')
   print "R0 [kpc]:", r0_kpc, x0_kpc-r0_kpc, x1_kpc-r0_kpc, '\n'
   # sigma
   sorted_sigma = sort(c.sigma)
   x = conf_inter(sorted_sigma, bestpars[4], 0.75, p)
   print "sigma:", bestpars[4], (x-bestpars[4]), '\n'
   # beta
   sorted_beta = sort(c.beta)
   x = conf_inter(sorted_beta, bestpars[5], 0.18, p)
   print "beta:", bestpars[5], (x-bestpars[5]), '\n'
