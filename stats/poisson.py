#!/usr/bin/env python

from numpy import *


def factorial(k):
   # calculate the factorial k! for *positive integers* or zero.
   if k < 0:
      raise ValueError, "Input needs to be positive!"
   if type(k) != type(0):
      raise ValueError, "Input is not a positive interger or zero!"
   if k == 0:
      return 1.0
   else:
      kf = 1
      while k > 0:
         kf = kf * k
         k = k - 1
      return float(kf)


def poisson_pdf(k, lamb):
   # calculate the Poisson distribution for the observed count k (integer) and
   # average count lamb
   k = int(k)
   p = lamb**k * exp(-1.*lamb) / factorial(k)
   return p


def poisson_cdf(k, lamb):
   # calculate the Poisson cumulative distribution function for the observed count k (integer)
   # and average count lamb
   k = int(k)
   cdf = 0.
   for ki in range(k+1):
      cdf = cdf + lamb**k / factorial(k)
   cdf = cdf * exp(-1.*lamb)
   return cdf


