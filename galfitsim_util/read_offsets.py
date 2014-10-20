#!/usr/bin/env python

import numpy as np
import os, sys, glob

def read_xyrange(offsets):
   f = open(offsets)
   lines = f.readlines()
   dx_array = []
   dy_array = []
   for i in range(len(lines)):
      if lines[i].startswith('#'):
         pass
      else:
         num, fname, offx, offy, xyrange, x, y, ncomp = lines[i].split()
         x_range, y_range = xyrange[1:-1].split(',')
         xmin, xmax = x_range.split(':')
         xmin = int(xmin); xmax = int(xmax)
         dx = xmax - xmin
         dx_array += [dx]
         ymin, ymax = y_range.split(':')
         ymin = int(ymin); ymax = int(ymax)
         dy = ymax - ymin
         dy_array += [dy]
   return np.array(dx_array), np.array(dy_array)

def read_xyrange_all(dir='.'):
   offsets = glob.glob(dir+'/offsets*')
   dx_array = np.zeros(0)
   dy_array = np.zeros(0)
   for off in offsets:
      dxr, dyr = read_xyrange(off)
      dx_array = np.concatenate([dx_array, dxr])
      dy_array = np.concatenate([dy_array, dyr])

   return dx_array, dy_array