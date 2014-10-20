#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from pygoods import Ftable, sextractor

class PlotCatalog(object):
   """
   A base class for making plots given a catalog. 
   Assumes the following attributes/methods:
   - self.c: this is the catalog instance
   (More...)
   """
   def get_col(self, colname, condition=None):
      assert hasattr(self.c, colname), "Column %s not in catalog..." % colname
      col = getattr(self.c, colname)
      if condition != None:
         return col[condition==True]
      else:
         return col

   def scatter_xy_cols(self, xcolname, ycolname, condition=None, ax=None, **scatterkw):
      """
      Makes a scatter plot give two column names.
      xcolname, ycolname: column names for each axis
      condition: criteria for selecting points to plot
      ax: a pyplot axis instance, default to None (create a new one)
      **scatterkw: additional keyword arguments fed directly to plt.scatter
      """
      xcol = self.get_col(xcolname, condition=condition)
      ycol = self.get_col(ycolname, condition=condition)
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      ax.scatter(xcol, ycol, **scatterkw)
      return ax

