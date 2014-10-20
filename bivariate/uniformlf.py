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

limits=array([[21.,27.],[-2.0,0.56]])
nbins=array([300,128]) 

def construct_model(r):
    dx = transfer.getdx(limits,nbins)
#   print dx
    y = int((r-limits[1,0])/dx[1])
    model = numarray.zeros(nbins,Float)
    model[:,y] = ones(shape(model[:,y]),Float)
    return model

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Usage: kernelfile.p r npts "
        sys.exit()

    kfile = sys.argv[1]
    f = open(kfile)
    kgrid = cPickle.load(f)

    r = string.atof(sys.argv[2])
    logr = log10(r)
    npts = string.atoi(sys.argv[3])

    model = construct_model(logr)
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


