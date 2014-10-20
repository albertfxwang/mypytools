#!/usr/bin/env python
# 
# Apply transfer function
#
# 7/11/04  v0.1
# v1.1.0 -- uses convolve.convolve2d instead of ndimage.correlate for speed
# 4/15/09 
# Modified by K-H Huang for using numpy and other improvements... 

from numpy import *
import pyfits
import sys
#from adaptivetransfer2d import *
#from scipy import signal
from scipy.ndimage import filters
from uniparam import limits,nbins,pixdx,dxgrid
from mlutil import findindex,bivmodel
import matplotlib.pyplot as plt
from mypytools import hconvolve
from multiprocessing import Process, Queue, Pool

#print dir(convolve)

__version__ = '1.1.0 4 January 2004' 
__author__ = 'Henry C. Ferguson, STScI'

class worker(object):
   def __init__(self, kgrid, bmodel):
      self.kgrid = kgrid
      self.bmodel = bmodel
   def __call__(self, key):
      return conv_kernel(key, self.kgrid, self.bmodel)
      

def assign_kernels(limits,nbins,kgrid,modpixdx=array([0.02,0.02])):
    """ Create an array of kernel indices (1-n) and a list of
        pointers to the corresponding kernels. Convert the
        kernels to numarrays in the process. 
        If modgrid is input, it returns the corresponding model
        grid that's consistent with kgrid.
        
        Returns:
	    kernel_list -- a list of numarray formatted kernels
            kernel_mask -- an array of indices into this list
                    This array has the dimensions of the *model* array.
                    A negative index means there is no kernel in this
                    portion of parameter space.
            model_list -- is empty if no input modgrid
    """
    ndim = len(nbins)
    xi = indices(nbins.tolist())*1.0
    xc = zeros(shape(xi)).astype('float') 
    kernelids = zeros((nbins),"int") 
    km = kernelids.ravel() # change of km will reference back to kernelids!!
    nkern = len(kernelids.ravel()) # nkern = total number of kernels in the grid
    dxgrid = kgrid.dxgrid
    die = 0
    for i in range(ndim):
        xc[i] = (xi[i]+0.5)*dxgrid[i]+limits[i,0] # centers of kernels
        ratgridsize = abs(modpixdx[i]-kgrid.pixdx[i])/modpixdx[i]
        if ratgridsize > 0.001:  # Error if kernel width differs by more than 0.1%
            print "Dimension %d: Kernel grid and model grid sizes disagree!" % (i)
            print "Kernel grid dx = %g ; Model grid dx = %g" % (kgrid.pixdx[i],modpixdx[i])
            die = 1
    if die:
        sys.exit(1)
#   print xc[0]
    kernel_list = []
    kernel_list_na = []
    for i in range(nkern):
        l = []
        for d in range(ndim):
            l = l+[xc[d].ravel()[i]]  
            # l[d] is the center coord of ith kernel in dth dimension
        a = array(l) 
        #print "getting kernel: ",a
	#print "i: ",i
	#print "xc: ",xc
        k = kgrid.getkernel(a)
        this_kernel_index = -1
        for kindex in range(len(kernel_list)):
            if k == kernel_list[kindex]:
                # if kernel k corresponding to center coord a is already in
                # kernel_list:
                this_kernel_index = kindex
                break
        if this_kernel_index < 0: # if k is not in kernel_list yet:
            # Convert kernel to numpy array and add to list
            kernel_list = kernel_list + [k]
            if k.completeness < 0:
                k.kernel = k.kernel*0.
            kna = array(k.kernel)
            kernel_list_na = kernel_list_na + [kna]
            this_kernel_index = len(kernel_list)-1
            if k.completeness < 0:
              print this_kernel_index,k.x0_bin, k.x1_bin, " kernel completeness: ",k.completeness, " setting to 0" 
            else:
              print this_kernel_index,k.x0_bin, k.x1_bin, " kernel completeness: ",k.completeness
	    #raw_input() 
        km[i] = this_kernel_index # modifies kernelids too!!
    kgrid.kernel_list = kernel_list_na  # contains all kernels in numpy arrays 
    kgrid.kernel_mask = kernelids
    return kernel_list_na,kernelids

def apply_transfer(bmodel,kgrid,min_comp=0.,
   fitsbasename="",fitshandle=0,fft=1):
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
   namodel = bmodel.model
   convolved_model = zeros(shape(namodel),"float")
   # Figure out which section in model to be convolved with which kernel
   nkey = 0
   for key in kgrid.kernels.keys():
      kern = kgrid.kernels[key]  # the kernel under consideration
      x0,x1 = kern.x0_cell,kern.x1_cell # the limits this kernel applies
      index0 = around((x0-limits[:,0])/pixdx).astype('int')
      index1 = around((x1-limits[:,0])/pixdx).astype('int')
      mid0,rid0 = index0
      mid1,rid1 = index1
      #print index0,index1
      if ((mid0<0) & (mid1<0)) | ((rid0<0) & (rid1<0)):
         # kernel outside model space 
         #print key, "kernel outside model limits"
         continue
      nkey += 1
      # if kernel bounds partially overlaps with model space: set indexes accordingly
      # so that the proper region in model is convolved with kernel
      if ((mid0<0) & (mid1>0)): mid0=0
      if ((rid0<0) & (rid1>0)): rid0=0
      if ((mid0>0) & (mid1>modshape[0])): mid1=modshape[0]
      if ((rid0>0) & (rid1>modshape[1])): rid1=modshape[1]
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
      # min_comp is used to mask out low completeness kernels
      if (sum(kern.kernel.ravel())>0) & (sum(section.ravel())>0):
         #if convmethod == signal.fftconvolve:
         #   convolved_sect = convmethod(section,kern.kernel,mode='full')
         #   convolved_sect = convolved_sect[kx0:-kx1,ky0:-ky1]
         #   # fftconvolve MUCH faster than convolve!!
         #elif convmethod == convolve.convolve2d:
         #   convolved_sect = convmethod(section,kern.kernel,mode='constant',fft=fft)
         #   #convolved_sect = convolved_sect[kx0:-kx1,ky0:-ky1]
         #else:
         convolved_sect = hconvolve.hconvolve(section,kern.kernel) # for hconvolve
      else:
         convolved_sect = 0.*section
      if fitshandle != 0:
         fitsname = "%s_%d_%d.fits" % (fitsbasename,key[0],key[1])
         print fitsname
         #fitshandle[0].data = convolved_sect[::-1,:]
         fitshandle.data = convolved_sect
         fitshandle.writeto(fitsname)
      convolved_model = convolved_model + convolved_sect
   bmodel = bivmodel(limits,pixdx,convolved_model)
   print "nkey", nkey
   return bmodel
 
def conv_kernel(key, kgrid, bmodel):
   # the worker function to perform convolution for a given kernel
   namodel = bmodel.model
   modshape = shape(namodel)
   limits = bmodel.limits
   pixdx = bmodel.pixdx
   kern = kgrid.kernels[key]  # the kernel under consideration
   x0,x1 = kern.x0_bin,kern.x1_bin # the limits this kernel applies
   index0 = around((x0-limits[:,0])/pixdx).astype('int')
   index1 = around((x1-limits[:,0])/pixdx).astype('int')
   mid0,rid0 = index0
   mid1,rid1 = index1
   #print index0,index1
   if ((mid0<0) & (mid1<0)) | ((rid0<0) & (rid1<0)):
      # kernel outside model space 
      #print key, "kernel outside model limits"
      return zeros(modshape)
   if ((mid0>modshape[0])&(mid1>modshape[0])) | ((rid0>modshape[1])&(rid1>modshape[1])):
      return zeros(modshape)
   # if kernel bounds partially overlaps with model space: set indexes accordingly
   # so that the proper region in model is convolved with kernel
   if ((mid0<0) & (mid1>0)): mid0=0
   if ((rid0<0) & (rid1>0)): rid0=0
   if ((mid0>0) & (mid1>modshape[0])): mid1=modshape[0]
   if ((rid0>0) & (rid1>modshape[1])): rid1=modshape[1]
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
   # min_comp is used to mask out low completeness kernels
   if (sum(kern.kernel.ravel())>0) & (sum(section.ravel())>0):
      convolved_sect = hconvolve.hconvolve(section,kern.kernel) # for hconvolve
   else:
      convolved_sect = 0.*section
   # put the results in the queues
   #qw_model.put(convolved_sect)
   #convolved_dict[key] = convolved_sect
   return convolved_sect

def apply_transfer_parallel(bmodel,kgrid,min_comp=0.,
   fitsbasename="",fitshandle=0,fft=1,convmethod=hconvolve.hconvolve,nproc=4):
   """ Apply transfer function (output of assign_kernels)
       Model must be a bivariate_lf.bivmodel object.
       Optionally output fits images of model convolved
       with the individual kernels.
       fitsbasename = root name of fits files
       fitshandle = template fits file handle already open
       Use python.multiprocessing module to make stepping through
       the kernels parallel.
   """
   modshape = shape(bmodel.model) 
   #print "nproc=", nproc
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
   namodel = bmodel.model
   convolved_model = zeros(shape(namodel),"float")
   #q_index = Queue()
   #q_model = Queue()

   # create the worker pool and execute the calculation
   pool = Pool()
   nsect = 0
   convolved_dict = {}
   
   # create the worker class
   wk = worker(kgrid, bmodel)
      
   p = pool.map_async(wk, kgrid.kernels.keys())   # pool.map_async is faster than pool.map
   results = p.get()
   #pool.close()
   #pool.join()   # join the pool AFTER emptying the queue
   #print "shape(results)", shape(results)

   convolved_model = reduce(lambda x,y:x+y, results)
   #convolved_model = array(results).sum(axis=0)   # slower than reduce...
   
   conv_model = bivmodel(limits,pixdx,convolved_model)
   return conv_model
 