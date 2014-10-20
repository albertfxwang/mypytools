#!/usr/bin/env python
#
# Version 1.0.0 H. Ferguson 9/18/03

"""Schechter function tools.
"""

__author__ =  'Henry C. Ferguson'
__version__=  '1.0.0'


import sys
from numpy import *
import random

def mag(f):
  return -2.5*log10(f)

def flux(m):
    return 10.**(-0.4*m)

def schechter(l,alpha):
    """ Return value of the Schechter function. L is in units L/L* """
    return exp(-l)*(l**alpha)

def mag_schechter(m,mstar,alpha):
    """ 
    Return value of the Schechter function (number per unit mag)
    for a given magnitude.
    """
    l = flux(m-mstar)
    try:
        phi = 0.4 * log(10.) * l**(alpha+1.) * exp(-l)
    except:
        print "Schecther overflow:"
        print "mstar = ",mstar
        print "alpha = ",alpha
        sys.exit() 
    return phi

def random_schechter(ntot,mstar,alpha,mfaint,mbright):
    """ Draw n random galaxies from a Schechter function. 
        Uses rejection method (not very efficient)"""
    lfaint = flux(mfaint-mstar)
    lbright = flux(mbright-mstar)
    if alpha > 0:
        print "random_schechter: alpha must be negative"
        sys.exit()
    smax = schechter(lfaint,alpha)
    result = zeros((ntot),"float")
    n = 0
    while n < ntot:
       x = random.uniform(lfaint,lbright)
       s = schechter(x,alpha)/smax
       y = random.random()    
       if y < s:
            result[n] = mag(x)+mstar
            n = n+1
    return result

if __name__ == "__main__":

# Test mag_schechter
#   mstar = float(sys.argv[1])
#   alpha = float(sys.argv[2])
#   m = arange(-24.,-10.,0.1)
#   s = mag_schechter(m,mstar,alpha)
#   for i in range(len(s)):
#        print m[i],s[i]

# Test random_schechter
    if len(sys.argv) < 6:
        print "usage schechter.py ntot mstar alpha mfaint mbright"
        sys.exit()
    ntot = int(sys.argv[1])
    mstar = float(sys.argv[2])
    alpha = float(sys.argv[3])
    mfaint = float(sys.argv[4])
    mbright = float(sys.argv[5])
    m = random_schechter(ntot,mstar,alpha,mfaint,mbright)
    for i in range(len(m)):
        print "%8.3f" % m[i]
