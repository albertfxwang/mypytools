#!/usr/bin/env python

from numpy import *
from pygoods import *


# Count the number of dropouts selected and eliminated by various
# criteria
# Visual inspection flags:
#  0 - normal-looking (isolated, not apparently disturbed, well detected, etc.)
#  1 - has stellar-like morphology or might not be resolved (will include in sample)
#  2 - confirmed a star
#  3 - multiple components (clumpy and/or merging) not deblended by SExtractor, therefore fitted together by GALFIT, OR same object broken up by SExtractor (will include in sample)
#  4 - confirmed a QSO
#  10 - diffraction spikes OR right on top of a diffraction spike
#  12 - near bright neighbor that photometry is unreliable
#  13 - not a galaxy but detected as a separate object by SExtractor for some reason 
#       OR lying too close to the edge to have reliable photometry OR saturated
#  14 - extremely LSB objects barely picked up by SExtractor (25 B-drops, 2.3%; 10 V-drops, 2.9%)
#  15 - with phot-z and/or spec-z below 2.0 so likely low-z interlopers (we only have phot-z information for GOODS-S sample and spec-z for very limited number of galaxies)
#  16 - this source is also covered by UDF. Since UDF is deeper, I will use UDF result.
#  17 - saturated star
#  18 - match with X-ray source; possible AGN
#  19 - very likely low-z object, but don't have a criterion to throw away
#  20 - source also exists in UDF, therefore use the UDF value

def count_bdrops(c1, c2, zlo=3.0):
   crit1 = (c1.cullflag != 20)
   crit2 = ones(len(c2), 'int')
   n1 = sum(crit1); n2 = sum(crit2)
   # step 1: count ALL sources selected by color
   print "Number selected by color from GOODS & HUDF are %d (%d in GOODS, %d in HUDF)" % (
      (n1+n2), n1, n2)
   # step 2: count ALL sources that pass visual inspection
   crit1 = crit1 & ((c1.cullflag==0)|(c1.cullflag==1)|(c1.cullflag==3)|(c1.cullflag==14)|\
      (c1.cullflag==19))
   crit2 = crit2 & ((c2.cullflag==0)|(c2.cullflag==1)|(c2.cullflag==3)|(c2.cullflag==14)|\
      (c2.cullflag==19))
   n1 = sum(crit1); n2 = sum(crit2)
   print "Number pass visual inspection is %d (%d in GOODS, %d in HUDF)" % ((n1+n2), n1, n2)
   # step 3: throw away sources with spec-z < zlo
   crit1 = crit1 & ((c1.specz<0.)|(c1.specz>=zlo))
   crit2 = crit2 & ((c2.specz<0.)|(c2.specz>=zlo))
   n1 = sum(crit1); n2 = sum(crit2)
   print "Number that either have no spec-z or have spec-z >= %.1f is %d (%d in GOODS, %d in HUDF)" % (zlo, (n1+n2), n1, n2)
   # step 4: throw away sources with 95% upper limit of phot-z <= zlo
   crit1 = crit1 & ((c1.photz_95max<0.)|(c1.photz_95max>=zlo))
   crit2 = crit2 & ((c2.photz_95max<0.)|(c2.photz_95max>=zlo))
   n1 = sum(crit1); n2 = sum(crit2)
   print "Number that either have no phot-z or have 95%% upper limit of phot-z >= %.1f is %d (%d in GOODS, %d in HUDF)" % (zlo, (n1+n2), n1, n2)
   return 0


def count_vdrops(c1, c2, zlo=4.0):
   crit1 = (c1.cullflag != 20)
   crit2 = ones(len(c2), 'int')
   n1 = sum(crit1); n2 = sum(crit2)
   # step 1: count ALL sources selected by color
   print "Number selected by color from GOODS & HUDF are %d (%d in GOODS, %d in HUDF)" % (
      (n1+n2), n1, n2)
   # step 2: count ALL sources that pass visual inspection
   crit1 = crit1 & ((c1.cullflag==0)|(c1.cullflag==1)|(c1.cullflag==3)|(c1.cullflag==14)|\
      (c1.cullflag==19))
   crit2 = crit2 & ((c2.cullflag==0)|(c2.cullflag==1)|(c2.cullflag==3)|(c2.cullflag==14)|\
      (c2.cullflag==19))
   n1 = sum(crit1); n2 = sum(crit2)
   print "Number pass visual inspection is %d (%d in GOODS, %d in HUDF)" % ((n1+n2), n1, n2)
   # step 3: throw away sources with spec-z < zlo
   crit1 = crit1 & ((c1.specz<0.)|(c1.specz>=zlo))
   crit2 = crit2 & ((c2.specz<0.)|(c2.specz>=zlo))
   n1 = sum(crit1); n2 = sum(crit2)
   print "Number that either have no spec-z or have spec-z >= %.1f is %d (%d in GOODS, %d in HUDF)" % (zlo, (n1+n2), n1, n2)
   # step 4: throw away sources with 95% upper limit of phot-z <= zlo
   crit1 = crit1 & ((c1.photz_95max<0.)|(c1.photz_95max>=zlo))
   crit2 = crit2 & ((c2.photz_95max<0.)|(c2.photz_95max>=zlo))
   n1 = sum(crit1); n2 = sum(crit2)
   print "Number that either have no phot-z or have 95%% upper limit of phot-z >= %.1f is %d (%d in GOODS, %d in HUDF)" % (zlo, (n1+n2), n1, n2)
   return 0
