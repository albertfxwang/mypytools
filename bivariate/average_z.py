#!/usr/bin/env python

from numpy import *
from pygoods import *
import dropout_sel as ds


def bestz(c, cullflags=[0,1,3,14,19], zlo=3.0):
   # construct best-redshift array AND best-redshift-upper-limit array
   bestz_arr = ones(len(c))*-1.
   bestzm_arr = ones(len(c))*-1.
   cullcrit = zeros(len(c),'int')
   for i in range(len(c)):
      if c.cullflag[i] in cullflags:
         cullcrit[i] = 1
      if c.specz[i] >= 0:
         bestz_arr[i] = c.specz[i]
         bestzm_arr[i] = c.specz[i]
      elif c.photz_weighted[i] > 0:
         bestz_arr[i] = c.photz_weighted[i]
         bestzm_arr[i] = c.photz_95max[i]
   bestz_arr = compress(cullcrit & (bestzm_arr >= zlo), bestz_arr)
   return bestz_arr

def average_median_z(c1, c2, cullflags=[0,1,3,14,19], zlo=3.0):
   bestz1 = bestz(c1, cullflags=cullflags, zlo=zlo)
   bestz2 = bestz(c2, cullflags=cullflags, zlo=zlo)
   bestz12 = concatenate((bestz1, bestz2))
   avg_z = average(bestz12)
   med_z = median(bestz12)
   return avg_z, med_z


def average_median_z_sim(c1, c2, drop='b'):
   if drop == 'b':
      dropcrit1 = ds.bdrops_sel(c1, zmagcut=26.5)[0]
      dropcrit2 = ds.bdrops_sel(c2, zmagcut=28.5)[0]
   elif drop == 'v':
      dropcrit1 = ds.vdrops_sel(c1, zmagcut=26.5, bstonlim=5.0)[0]
      dropcrit2 = ds.vdrops_sel(c2, zmagcut=28.5, bstonlim=5.0)[0]
   bestz1 = compress(dropcrit1, c1.z_input)
   bestz2 = compress(dropcrit2, c2.z_input)
   bestz = concatenate((bestz1, bestz2))
   avg_z = average(bestz)
   med_z = median(bestz)
   return avg_z, med_z


