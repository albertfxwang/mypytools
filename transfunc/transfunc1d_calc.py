#!/usr/bin/env python

from numpy import *
from pygoods import *
import KPDFadaptnumpy as KPDF
import draw_ndist as dn
import matplotlib as mpl
import matplotlib.pyplot as plt
import cPickle
import mlutil
from scipy.interpolate import RectBivariateSpline

"""
A general-purpose script that calculates 1D transfer-function kernels in cells 
of the parameter space, using data from simulation output.
"""

class kernel1d(object):
   """ A convolution kernel associated with a certain range of parameter space
   """
   def __init__(self,x0_cell,x1_cell,pixdx,completeness,kernel_array):
      self.x0_cell = x0_cell                         
      # Where to start using this kernel
      self.x1_cell = x1_cell                         
      # Where to stop using this kernel
      self.xc_cell = (x0_cell + x1_cell) / 2.
      # The coordinate at the center of this cell
      self.pixdx = array(pixdx)                         
      # kernel pixel size (a length-2 array)
      self.i0 = -array(shape(kernel_array))/2    
      # the pixel index of the lower-left corner, relative to center
      self.i1 = array(shape(kernel_array))/2     
      # the pixel index of the upper-right corner, relative to center
      self.rx0 = -array(shape(kernel_array))*pixdx/2.  
      # the pixel coordinate of the lower-left corner, relative to center
      self.rx1 = array(shape(kernel_array))*pixdx/2.   
      # the pixel coordinate of the upper-right corner, relative to center
      self.completeness = completeness   
      # completeness is specified separately from kernel_array
      self.kernel = kernel_array   
      # always normalized to sum(kernel_array.ravel()) = 1

   def update(self,kernel_array,pixdx=None):
      if pixdx == None:
         pass
      else:
         self.pixdx = pixdx
      self.i0 = -array(shape(kernel_array))/2
      self.i1 = array(shape(kernel_array))/2
      self.rx0 = -array(shape(kernel_array))*self.pixdx/2.
      self.rx1 = array(shape(kernel_array))*self.pixdx/2.
      #self.completeness = sum(kernel_array.ravel()) * self.pixdx.prod()
      self.kernel = kernel_array

   