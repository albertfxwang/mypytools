#!/usr/bin/env python
""" Iterative rejection of an array with a ceiling and a floor """

__author__ = "Henry C. Ferguson, Space Telescope Science Institute"
__copyright__ = "Copyright 2012, Henry C. Ferguson"
__credits__ = ["Henry C. Ferguson"]
__license__ = "FreeBSD"
__version__ = "0.1"
__maintainer__ = "Henry C. Ferguson"
__email__ = "ferguson@stsci.edu"
__status__ = "Development"

from numpy import *

def iterstat(x,niter,low=None,high=None,ceiling=None,floor=None,justmean=False,verbose=False):
    """ Statistics with iterative rejection and a ceiling and a floor

        Arguments:
        x -- Array for which to determine the clipped mean & standard deviation
        niter -- Number of iterations
        low -- reject points lower than this number of standard deviations on each iteration
        high -- reject points high than this number of standard deviations on each iteration
        ceiling -- reject points higher than this value
        floor -- reject points higher than this value
    """

    if verbose:
        print "%4s %10s %12s %12s %12s" % ('iter','npix','median','mean','std')

    # Check if the inputs are reasonable
    if ceiling and floor:
       if ceiling <= floor:
          print "Iterstat: ceiling < floor; ignoring both of them"
	  ceiling=None
	  floor=None

    for n in range(niter):
        m = median(x)
        s = x.std()
        filt = ones(x.shape,dtype=bool)
        if low != None:
            filt = filt & (x>m-low*s)
        if high != None:
            filt = filt & (x<m+high*s) 
        if ceiling != None:
            filt = filt & (x<ceiling)
        if floor != None:
            filt = filt & (x>floor)
        if filt.sum() < 1:  # If we are left with no points, just return previous value
	   print "Iterstat: All pixels would be rejected, returning with no rejection"
           break
        keep = where(filt)
        x = x[keep]
	if verbose:
            print "%4d %10d %12.3g %12.3g %12.3g" % (n,filt.sum(),median(x),x.mean(),x.std())
    if justmean:
        return m
    else:
        return m,s,len(x.flat)
