#!/usr/bin/env python

from numpy import *
from pygoods import *
import bivariate_lf as bl
import bivariate_fit as bf
import fit_lbg as fl
import mlutil

parb = array([-1.61711035, -20.53430348, 0.81505052, 0.76553555, 0.20435599])
parv = array([-1.68527018, -20.50967054, 0.8757289, 0.85187255, 0.26594964])
kgrid1 = mlutil.readkgrid('kernel_I.p')
kgrid2 = mlutil.readkgrid('kernel_I_udf.p')
kgrid3 = mlutil.readkgrid('kernel_Z.p')
kgrid4 = mlutil.readkgrid('kernel_Z_udf.p')
mci = bl.mconvert('M1500_to_i.txt')
mcz = bl.mconvert('M1500_to_z.txt')

def logl_sigma_bdrops(sarr=arange(0.3,1.0,0.05)):
   mag1, re1, crit1 = fl.cleandata('bdrops_gf_v2.cat',chisqnulim=0.4,magautolim=26.5,
      limits=bl.limits1,drop='b')
   mag2, re2, crit2 = fl.cleandata('bdrops_udf_gf_v2.cat',chisqnulim=5.0,magautolim=28.5,
      limits=bl.limits2,drop='b')
   data1 = array([mag1, log10(re1)])
   data2 = array([mag2, log10(re2)])
   loglarr = zeros(len(sarr))
   for i in range(len(sarr)):
      par = parb.copy()
      par[3] = sarr[i]
      loglarr[i] = bf.mlfunc(par,data1,data2,bl.limits1,bl.limits2,bl.pixdx,
         kgrid1,kgrid2,1.0,1.0,-21.0,'zdgrid_bdrops.p','zdgrid_bdrops_udf.p',
         'M1500_to_i.txt',1,mci(4.0),0,-1.,'phistar','b')
   return loglarr
   
def logl_sigma_vdrops(sarr=arange(0.3,1.0,0.05)):
   mag1, re1, crit1 = fl.cleandata('vdrops_gf_v2.cat',chisqnulim=0.5,magautolim=26.5,
      limits=bl.limits1,drop='v')
   mag2, re2, crit2 = fl.cleandata('vdrops_udf_gf_v2.cat',chisqnulim=5.0,magautolim=28.5,
      limits=bl.limits2,drop='v')
   data1 = array([mag1, log10(re1)])
   data2 = array([mag2, log10(re2)])
   loglarr = zeros(len(sarr))
   for i in range(len(sarr)):
      par = parv.copy()
      par[3] = sarr[i]
      loglarr[i] = bf.mlfunc(par,data1,data2,bl.limits1,bl.limits2,bl.pixdx,
         kgrid1,kgrid2,1.0,1.0,-21.0,'zdgrid_vdrops.p','zdgrid_vdrops_udf.p',
         'M1500_to_z.txt',1,mcz(5.0),0,-1.,'phistar','v')
   return loglarr