#!/usr/bin/env python

from numpy import *
from pygoods import *
import pyfits
import sys, os, time
import mlutil, bivmodel
import hconvolve
from multiprocessing import Process, Queue, Pool, Manager
from transfunc import transfunc_calc
import cPickle

__version__ = "1.0"
__author__ = "Kuang-Han Huang / H. C. Ferguson"

def readkgrid(kgridfile):
    f = open(kgridfile)
    kgrid = cPickle.load(f)
    f.close()
    return kgrid
   
def conv_kernel(key, kgrid, bmodel):
   """
   Given a kernel grid kgrid, find the corresponding kernel using the provided key.
   Then convolve the corresponding region of bmodel with the kernel and returned 
   the convolved model array.
   key -- a list of indices for a given kernel
   kgrid -- the transfunc_calc.kernelgrid instance
   bmodel -- a bivmodel.bivmodel instance
   """
   # the worker function to perform convolution for a given kernel
   namodel = bmodel.model
   modshape = shape(namodel); print modshape
   limits = bmodel.limits
   pixdx = bmodel.pixdx
   kern = kgrid.kernels[key]  # the kernel under consideration
   print shape(kern.kernel)
   x0,x1 = kern.x0_cell,kern.x1_cell # the limits this kernel applies
   index0 = around((x0-limits[:,0])/pixdx).astype('int')
   index1 = around((x1-limits[:,0])/pixdx).astype('int')
   mid0,rid0 = index0
   mid1,rid1 = index1
   #print "mid0,mid1,rid0,rid1",mid0,mid1,rid0,rid1
   # if this kernel is COMPLETELY outside of the limits -- do not perform convolution
   if (mid1<0) | (rid1<0):
      return zeros(modshape)
   if (mid0>modshape[0]) | (rid0>modshape[1]):
      return zeros(modshape)
   # if kernel bounds partially overlaps with model space: set indexes accordingly
   # so that the proper region in model is convolved with kernel
   mid0 = maximum(mid0,0)
   rid0 = maximum(rid0,0)
   mid1 = minimum(mid1,modshape[0])
   rid1 = minimum(rid1,modshape[1])
   #print "mid0,mid1,rid0,rid1",mid0,mid1,rid0,rid1
   kx,ky = shape(kern.kernel) # the kernel array dimension
   if kx%2 == 0:
      kx0 = kx/2 - 1    # padding before section
      kx1 = kx/2        # padding after section
   else:
      kx0 = kx/2
      kx1 = kx/2
   if ky%2 == 0:
      ky0 = ky/2 - 1
      ky1 = ky/2
   else:
      ky0 = ky/2
      ky1 = ky/2
   section = zeros(shape(namodel),"float")
   section[mid0:mid1,rid0:rid1] = namodel[mid0:mid1,rid0:rid1]  # paste the model onto empty image   
   # min_comp is used to mask out low completeness kernels --- not implemented now
   if (sum(kern.kernel.ravel())>0) & (sum(section.ravel())>0):
      convolved_sect = hconvolve.hconvolve(section,kern.kernel) # for hconvolve
   else:
      convolved_sect = 0.*section
   return convolved_sect

class worker(object):
   def __init__(self, kgrid, bmodel):
      self.kgrid = kgrid
      self.bmodel = bmodel
   def __call__(self, key):
      return conv_kernel(key, self.kgrid, self.bmodel)
 
def apply_transfer_2d(bmodel, kgrid, nproc_model=1):
   """ Apply transfer function (output of assign_kernels)
       Model must be a bivariate_lf.bivmodel object.
       Optionally output fits images of model convolved
       with the individual kernels.
       fitsbasename = root name of fits files
       fitshandle = template fits file handle already open
   """
   modshape = shape(bmodel.model) 
   dd = kgrid.pixdx/1.e3   # tolerance for index determination
   limits = bmodel.limits
   pixdx = bmodel.pixdx
   compdx = (abs(pixdx-kgrid.pixdx)<=dd)
   if not compdx.all():
      print "Model and kernel pixel size not compatible:"
      print "model dx = ", pixdx
      print "kernel dx = ", kgrid.pixdx
      sys.exit(1)
   # construct the convolved model grid
   # Multiply the entire model by GALFIT completeness before convolution
   # Then all kernels are normalized to sum to 1.0
   namodel = bmodel.model * kgrid.completeness_model
   convolved_model = zeros(shape(namodel),"float")
   # Figure out which section in model to be convolved with which kernel
   chunksize = len(kgrid.kernels.keys()) / nproc_model
   manager = Manager()
   q_out = manager.Queue()
   processes = []
   for i in range(nproc_model):
      i0 = i * chunksize
      i1 = minimum((i+1)*chunksize, len(kgrid.kernels.keys()))
      klist = kgrid.kernels.keys()[i0:i1]
      p = Process(target=apply_transfer_2d_worker, 
                     args=(bmodel, kgrid, klist, q_out))
      processes += [p]
      p.start()
   # join the processes
   for i in range(nproc_model):
      processes[i].join()
   # collect output
   while not q_out.empty():
      convolved_model = convolved_model + q_out.get()
   bmodel = bivmodel.bivmodel(limits,pixdx,convolved_model)
   return bmodel 

def apply_transfer_2d_worker(bmodel, kgrid, klist, q_out):
   for k in klist:
      kern = kgrid.kernels[k]  # the kernel under consideration
      #print shape(kern.kernel)
      x0,x1 = kern.x0_cell,kern.x1_cell # the limits this kernel applies
      index0 = around((x0-bmodel.limits[:,0])/bmodel.pixdx).astype('int')
      index1 = around((x1-bmodel.limits[:,0])/bmodel.pixdx).astype('int')
      mid0,rid0 = index0
      mid1,rid1 = index1
      modshape = shape(bmodel.model) 
      namodel = bmodel.model
      #print "mid0,mid1,rid0,rid1",mid0,mid1,rid0,rid1
      # if this kernel is COMPLETELY outside of the limits -- do not perform convolution
      if (mid1<0) | (rid1<0):
         continue
      elif (mid0>modshape[0]) | (rid0>modshape[1]):
         continue
      # if kernel bounds partially overlaps with model space: set indexes accordingly
      # so that the proper region in model is convolved with kernel
      mid0 = maximum(mid0,0)
      rid0 = maximum(rid0,0)
      mid1 = minimum(mid1,modshape[0])
      rid1 = minimum(rid1,modshape[1])
      kx,ky = shape(kern.kernel) # the kernel array dimension
      if kx%2 == 0:
         kx0 = kx/2 - 1    # padding before section
         kx1 = kx/2        # padding after section
      else:
         kx0 = kx/2
         kx1 = kx/2
      if ky%2 == 0:
         ky0 = ky/2 - 1
         ky1 = ky/2
      else:
         ky0 = ky/2
         ky1 = ky/2
      section = zeros(shape(namodel),"float")
      section[mid0:mid1,rid0:rid1] = namodel[mid0:mid1,rid0:rid1]  
      # paste the model onto empty image
      if (sum(kern.kernel.ravel())>0) & (sum(section.ravel())>0):
         convolved_sect = hconvolve.hconvolve(section,kern.kernel) # for hconvolve
      else:
         convolved_sect = 0.*section
      q_out.put(convolved_sect)
      #convolved_model = convolved_model + convolved_sect

   
def apply_transfer_2d_parallel(bmodel,kgrid,nproc=1):
   """ 
   Apply transfer function (output of assign_kernels)
   Model must be a bivariate_lf.bivmodel object.
   Optionally output fits images of model convolved
   with the individual kernels.
   fitsbasename = root name of fits files
   fitshandle = template fits file handle already open
   Use python.multiprocessing module to make stepping through
   the kernels parallel.
   """
   modshape = shape(bmodel.model) 
   dd = kgrid.pixdx/1.e3   # tolerance for index determination
   limits = bmodel.limits
   pixdx = bmodel.pixdx
   compdx = (abs(pixdx-kgrid.pixdx)<=dd)
   if not compdx.all():
      print "Model and kernel pixel size not compatible:"
      print "model dx = ", pixdx
      print "kernel dx = ", kgrid.pixdx
      sys.exit(1)
      
   # construct the convolved model grid
   convolved_model = zeros(shape(bmodel.model),"float")
   #q_index = Queue()
   #q_model = Queue()

   # create the worker pool and execute the calculation
   pool = Pool()
   nsect = 0
   convolved_dict = {}
   
   # create the worker class
   wk = worker(kgrid, bmodel)
   print "start convolving..."
   p = pool.map_async(wk, kgrid.kernels.keys())   # pool.map_async is faster than pool.map
   results = p.get()
   #pool.close()
   #pool.join()   # join the pool AFTER emptying the queue
   #print "shape(results)", shape(results)
   print "finished convolving..."
   convolved_model = reduce(lambda x,y:x+y, results)
   
   conv_model = bivmodel.bivmodel(limits,pixdx,convolved_model)
   return conv_model
