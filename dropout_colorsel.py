#!/usr/bin/env python

from numpy import *
from pygoods import *


# the collection of dropout-selection functions from the literature

goods_acszpt = {'B':25.65288, 'V':26.49341, 'I':25.64053, 'Z':24.84315}
udf09_wfc3zpt = {'Y':26.261, 'J':26.237, 'H':25.932}


def flux2mag(flux,zpt):
   # calculate magnitude given flux and zeropoint
   mag = zpt - 2.5 * log10(flux)
   return mag


def flux_uplim(fluxarr, fluxerrarr, zpt, magarr):
   # replace the flux with upper limits when S/N < 1
   for i in range(len(fluxarr)):
      if fluxarr[i] / fluxerrarr[i] <= 1.0:
         if fluxarr[i] > 0.:
            magarr[i] = flux2mag(fluxerrarr[i], zpt)
         else:
            magarr[i] = 99.0
   return magarr


def bdrops_colorsel_g04(bmags,vmags,zmags):
   """Giavalisco et al. 2004 criteria"""
   # make sure to use upper limits for S/N < 1 in b-band and v-band
   bmv = bmags - vmags
   vmz = vmags - zmags
   c1 = (bmv >= (1.2 + 1.4 * vmz))
   c2 = (bmv >= 1.2)
   c3 = (vmz <= 1.2)
   return (c1 & c2 & c3)


def vdrops_colorsel_g04(vmags,imags,zmags):
   """Giavalisco et al. 2004 V-dropouts color selection
      colorsel_g04(vmags,imags,zmags)
   """
   vmi = vmags - imags
   imz = imags - zmags
   c1 = (vmi > 1.5 + 0.9 * imz)
   c2 = (vmi > 2.0)
   c3 = (vmi >= 1.2)
   c4 = (imz <= 1.3)
   return  ((c1 | c2) & c3 & c4)


def ydrops_colorsel_b11(ymags, jmags, hmags, ston_dic, opt_snlim=2.0):
   """Bouwens et al. 2011, ApJ, 737, 90 Y-dropout (z~8) selection criteria
      Requires WFC3/IR colors as well as signal-to-noise ratios in the optical bands (BViz).
   """
   ymj = ymags - jmags
   jmh = jmags - hmags
   c1 = (ymj > 0.45)
   c2 = (jmh < 0.5)
   #y_ston = (ston_dic['f105w'] > 3.0)
   j_ston = (ston_dic['J'] > 3.5)
   h_ston = (ston_dic['H'] > 3.0)
   ir_stoncrit = (h_ston & j_ston)
   opt_stoncrit1 = ((ston_dic['B']<opt_snlim) & (ston_dic['V']<opt_snlim) &\
      (ston_dic['I']<opt_snlim) & (ston_dic['Z']<opt_snlim))
   opt_stoncrit2 = (((ston_dic['B']>1.5)+(ston_dic['V']>1.5)+(ston_dic['I']>1.5)+\
      (ston_dic['Z']>1.5)) <= 1)
   return (c1 & c2 & ir_stoncrit & opt_stoncrit1 & opt_stoncrit2)


def ydrops_colorsel_yan11(ymags, jmags, hmags, jston_auto, jston_iso, 
   hston_auto, hston_iso, ston_dic, opt_snlim=2.0):
   # the Y-dropout criteria from Haojing Yan
   # as used in Yan et al. 2012
   
   ymj = ymags - jmags
   jmh = jmags - hmags
   c1 = (ymj >= 0.8)
   c2 = (jmh <= 0.3)
   j_ston = (jston_auto >= 5.0) | (jston_iso >= 5.0)
   h_ston = (hston_auto >= 5.0) | (hston_iso >= 5.0)
   ir_stoncrit = (j_ston & h_ston)
   opt_stoncrit1 = ((ston_dic['b']<opt_snlim) & (ston_dic['v']<opt_snlim) &\
      (ston_dic['i']<opt_snlim) & (ston_dic['z']<opt_snlim))
   return (c1 & c2 & ir_stoncrit & opt_stoncrit1)


