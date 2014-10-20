#!/usr/bin/env python
#
# Construct a simulated data set from a transfer function.
# 
# Inputs: kernelfile.p magnitude re npoints
# Output: magnitudes, re  for n galaxies drawn from the appropriate
#      transfer function kernel
#

import cPickle
import transfer
from apply_transfer_nd import *
from mlutil import *
from Numeric import *
from numprint import *
import numarray
import pyfits
import glob, os, sys, string

__author__ =  'Henry C. Ferguson'
__version__=  '0.1.0'

limits=array([[21.,27.],[-1.3,0.7]])
nbins=array([300,200]) 

def construct_model(mag,r):
    dx = transfer.getdx(limits,nbins)
#   print dx
    x = int((mag-limits[0,0])/dx[0])
    y = int((r-limits[1,0])/dx[1])
    model = numarray.zeros(nbins,Float)
    model[x,y] = 1.0
    return model

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print "Usage: kernelfile.p mag r npts "
        sys.exit()

    kfile = sys.argv[1]
    f = open(kfile)
    kgrid = cPickle.load(f)

    mag = string.atof(sys.argv[2])
    r = string.atof(sys.argv[3])
    logr = log10(r)
    npts = string.atoi(sys.argv[4])

    model = construct_model(mag,logr)
    kernel_list, kernel_mask = assign_kernels(limits,nbins,kgrid)
#   print shape(kernel_list)
#   print shape(kernel_mask)
#   print kernel_mask[5,5]
#   print kernel_mask[15,25]
#   for i in range(len(kernel_list)):
#       print "shape(kernel_list[%d]):" % i,
#       print shape(kernel_list[i])
    smeared_model = apply_transfer(model,kernel_mask,kernel_list)
#   print "shape(smeared_model) ",shape(smeared_model)
#   print "limits ",limits
#   print "nbins ",nbins
    nsmeared_model = array(smeared_model.tolist())
    (magout,y) = draw_from_pdf(npts,nsmeared_model,limits,nbins)
#   print magout
#   print y
    ten = 10.*ones(shape(y))
    radius = ten**y

    l = format("%8.3f %8.3f",magout,radius)
    print l


