#!/usr/bin/env python

from numpy import *
import pysynphot as S
from pysynphot import Extinction
import matplotlib as mpl
import matplotlib.pyplot as plt
from sedtools import magredshift as mr
import os

biv = os.getenv('biv')
sp = S.FileSpectrum(biv+'/bivariate_fit/simcatalogs/idrops/lbgtemp_z6_0.4Zsolar_constSFH_100Myr.sed')
sp_ext = sp * Extinction(0.15, 'xgal')
#zarr = arange(5.5, 7.0, 0.4)
zarr = [5.5, 6.0, 6.5, 7.0]
f435w = S.ObsBandpass('acs,wfc2,f435w')
f606w = S.ObsBandpass('acs,wfc2,f606w')
f775w = S.ObsBandpass('acs,wfc2,f775w')
f850lp = S.ObsBandpass('acs,wfc2,f850lp')
f125w = S.ObsBandpass('wfc3,ir,f125w')
f160w = S.ObsBandpass('wfc3,ir,f160w')
filters = [f435w, f606w, f775w, f850lp, f125w, f160w]

fig = plt.figure()
plt.yscale('log')

for f in filters:
   plt.plot(f.wave, f.throughput/sum(f.throughput), color='0.3')

for z in zarr:
   sp_red = mr.mag_redshift(sp_ext, z, f160w)[1]
   plt.plot(sp_red.wave, 50.*sp_red.flux/sum(sp_red.flux), label='z=%.1f'%z, 
            lw=1.5)

plt.legend(loc=2)
plt.ylim(1.e-5, 1.e-1)
plt.xlim(1000., 20000.)
