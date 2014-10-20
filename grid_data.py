#!/usr/bin/env python

import numpy as np

"""
Grid data points into regular bins. Now implementing both the 1D and 2D cases.
Note that the array for bins should contain the upper (right) edge of the 
last bin.
"""

def grid_data_1d(x, binedges, option=1):
   """
   Return the elements in x that satisfy (x>=binedge[0]) and (x<binedge[1]).
   Should be as fast as possible for a large array.
   """
   # Method 1: directo slicing using array index
   if option==1:
      y = x[(x>=binedges[0]) & (x<binedges[1])]
   # Method 2: use searchsorted
   elif option==2:
      z = np.searchsorted(binedges, x)
      y = x[z==1]
   return y


def grid_data_2d(x, y, xbinedges, ybinedges, option=1):
   """
   Return an array of shape (2, N), where N is the number of elements in x
   that are within xbinedges and in y that are within ybinedges.
   """
   if option==1:
      crit = [(x>=xbinedges[0])&(x<xbinedges[1])&(y>=ybinedges[0])&(y<ybinedges[1])]
      z = np.array([x[crit], y[crit]])
   elif option==2:
      x1 = np.searchsorted(xbinedges, x)
      y1 = np.searchsorted(ybinedges, y)
      z1 = ((x1==1) & (y1==1))
      z = np.array([x[z1], y[z1]])
   elif option==3:
      crit = (x>=xbinedges[0])&(x<xbinedges[1])
      crit = crit & (y>=ybinedges[0])&(y<ybinedges[1])
      z = np.array([x[crit], y[crit]])
   return z