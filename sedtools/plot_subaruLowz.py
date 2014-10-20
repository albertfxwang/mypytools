#!/usr/bin/env python

import pylab as plt
import galtemp_colors as gc
import pysynphot as S
import numpy as np
from pygoods import sextractor

tempdir = '/Users/khuang/Dropbox/codes/mypytools/sedtools/galtemplates'
## Coleman, Wu & Weedman (1980) local galaxy templates
Sbc = tempdir + '/Gsbc.spec'
Scd = tempdir + '/Gscd.spec'
ES0 = tempdir + '/Geso.spec'
Imm = tempdir + '/Gimm.spec'
## LBG template with const. SFH & age = 200 Myr, solar Z, Salpeter IMF
LBG = tempdir + '/lbgtemp_z3_Zsolar_constSFH_200Myr.sed'

localTemplates = {'Sbc':Sbc, 'Scd':Scd, 'ES0':ES0, 'Imm':Imm, 'LBG':LBG}

## Subaru filters
bpdir = '/Users/khuang/Dropbox/codes/mypytools/sedtools/bandpass/CANDELS_FILTERS'
iSubaru = S.FileBandpass(bpdir + '/SuprimeCAM/i_subaru.dat')
iSubaru.name = 'iSubaru'
jSubaru = S.FileBandpass(bpdir + '/MOIRCS/J_MOIRCS.dat')
jSubaru.name = 'jSubaru'
hSubaru = S.FileBandpass(bpdir + '/MOIRCS/H_MOIRCS.dat')
hSubaru.name = 'hSubaru'
kSubaru = S.FileBandpass(bpdir + '/MOIRCS/K_MOIRCS.dat')
kSubaru.name = 'kSubaru'
subFilters = [iSubaru, jSubaru, hSubaru, kSubaru]

def calcTempMagnitudes(templates=localTemplates, zarray=np.arange(0.02,3.02,0.02)):
   for spec in templates:
      filename = '%s_subaruiJHK.mag' % spec
      print "Working on %s..." % filename
      g = gc.GalTemplateColors(localTemplates[spec], filters=subFilters)
      g.writeMagnitudes(zarray, magfile=filename)

def plotTempColors(templates=localTemplates):
   fig = plt.figure(figsize=(12,10))
   ax1 = fig.add_subplot(2,2,1)  # i - J
   ax2 = fig.add_subplot(2,2,2)  # i - H
   ax3 = fig.add_subplot(2,2,3)  # i - K
   for spec in templates:
      filename = '%s_subaruiJHK.mag' % spec
      c = sextractor(filename)
      ax1.plot(c.z, c.isubaru - c.jsubaru, label=spec)
      ax2.plot(c.z, c.isubaru - c.hsubaru, label=spec)
      ax3.plot(c.z, c.isubaru - c.ksubaru, label=spec)

   for ax in [ax1, ax2, ax3]:
      ax.set_xlabel('Redshift')
      ax.legend(loc=2)
      ax.set_ylim(-1.,6.)
   ax1.set_ylabel(r'$i - J$')
   ax2.set_ylabel(r'$i - H$')
   ax3.set_ylabel(r'$i - K$')

