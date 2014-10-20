#!/usr/bin/env python

import numpy as N
from pygoods import *
import pyfits
import scipy
import mlutil
import transfunc_apply, transfunc_calc
import matplotlib as mpl
import matplotlib.pyplot as plt
import cPickle

class bivmodel(object):
   """
   An object that contains information of model limits, pixel size, etc.
   """
   def __init__(self, limits, pixdx, model_array, axeslabel=['magnitude','log10(Re)']):
      self.limits = limits
      self.pixdx = pixdx
      self.model = model_array
      self.axeslabel = axeslabel
      mags = N.arange(limits[0,0],limits[0,1],pixdx[0])
      logr = N.arange(limits[1,0],limits[1,1],pixdx[1])
      coords = N.zeros((len(mags),len(logr),2))
      for i in range(len(mags)):
         coords[i,:,0] = mags[i]
      for j in range(len(logr)):
         coords[:,j,1] = logr[j]
      self.coords = coords  
      # self.coords[i,j] returns (mag,logr) corresponding to the (i,j)th
      # pixel
      
   def copy(self):
      x = self.__init__(self.limits,self.pixdx,self.model)
      return x
        
   def apply_transfer_2d(self,kgrid=None):
      """
      calls transfunc_apply.apply_transfer_2d().
      """
      if kgrid==None: kgrid = self.kgrid
      #conv_model = transfunc_apply.apply_transfer_2d_parallel(self,kgrid,nproc=nproc)
      conv_model = transfunc_apply.apply_transfer_2d(self,kgrid)
      self.__dict__.update(conv_model.__dict__)  # update itself
      
   def show(self,tick_widths=None,imshow_kw={}):
      if tick_widths == None: tick_widths = array([0.5,0.2])
      # figure out the coordinates of the ticks
      xtickcoords = arange(self.limits[0,0],self.limits[0,1],tick_widths[0])
      xticks = (xtickcoords-self.limits[0,0]) / self.pixdx[0]
      ytickcoords = arange(self.limits[1,0],self.limits[1,1],tick_widths[1])
      yticks = (ytickcoords-self.limits[1,0]) / self.pixdx[1]
      plt.imshow(self.model.swapaxes(0,1), origin='lower', **imshow_kw)
      plt.xticks(xticks,xtickcoords)
      plt.yticks(yticks,ytickcoords)
      plt.xlabel(self.axeslabel[0])
      plt.ylabel(self.axeslabel[1])
      
      
      
      