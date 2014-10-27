#!/usr/bin/env python

import numpy as np
import FileSED
import pysynphot as S
import matplotlib.pyplot as plt
from pygoods import sextractor
import cosmoclass
import kcorr
# For plotting colors of stars, look at starlib.py

cosmo_def = cosmoclass.cosmoclass(H0=70., omega_m=0.3, omega_l=0.7)
# lbgtemp = 'galtemplates/lbgtemp_z6_0.4Zsolar_constSFH_100Myr.sed'
# Below are HST ACS/WFC3IR filters in CLASH; likely all I'll need for high-z 
# stuff?
banddir = '/Users/khuang/Dropbox/codes/mypytools/sedtools/bandpass'
uvimos = S.FileBandpass('%s/U_vimos.bp' % banddir)
uvimos.name = 'UVIMOS'
f435w = S.ObsBandpass('acs,wfc2,f435w')
f475w = S.ObsBandpass('acs,wfc2,f475w')
f555w = S.ObsBandpass('acs,wfc2,f555w')
f606w = S.ObsBandpass('acs,wfc2,f606w')
f625w = S.ObsBandpass('acs,wfc2,f625w')
f775w = S.ObsBandpass('acs,wfc2,f775w')
f814w = S.ObsBandpass('acs,wfc2,f814w')
f850lp = S.ObsBandpass('acs,wfc2,f850lp')
f098m = S.ObsBandpass('wfc3,ir,f098m')
f105w = S.ObsBandpass('wfc3,ir,f105w')
f110w = S.ObsBandpass('wfc3,ir,f110w')
f125w = S.ObsBandpass('wfc3,ir,f125w')
f140w = S.ObsBandpass('wfc3,ir,f140w')
f160w = S.ObsBandpass('wfc3,ir,f160w')
VJohnson = S.FileBandpass('%s/CANDELS_FILTERS/Johnson/Johnson_V.dat' % banddir)
VJohnson.name = 'VJohnson'
IJohnson = S.FileBandpass('%s/CANDELS_FILTERS/Johnson/Ifilter.dat' % banddir)
IJohnson.name = 'IJohnson'
iSubaru = S.FileBandpass('%s/CANDELS_FILTERS/SuprimeCAM/i_subaru.dat' % banddir)
iSubaru.name = 'iSubaru'
jSubaru = S.FileBandpass('%s/CANDELS_FILTERS/MOIRCS/J_MOIRCS.dat' % banddir)
jSubaru.name = 'jSubaru'
hSubaru = S.FileBandpass('%s/CANDELS_FILTERS/MOIRCS/H_MOIRCS.dat' % banddir)
hSubaru.name = 'hSubaru'
kSubaru = S.FileBandpass('%s/CANDELS_FILTERS/MOIRCS/K_MOIRCS.dat' % banddir)
kSubaru.name = 'kSubaru'
IRAC1 = S.FileBandpass('%s/CANDELS_FILTERS/IRAC/irac_ch1.dat' % banddir)
IRAC2 = S.FileBandpass('%s/CANDELS_FILTERS/IRAC/irac_ch2.dat' % banddir)
# default_filters = [f435w, f814w, f098m, f125w, f160w]
default_filters = [uvimos, f435w, f475w, f555w, f606w, f625w, f775w, f814w, \
   f850lp, f098m, f105w, f110w, f125w, f140w, f160w]
# local galaxy templates
SB2 = 'galtemplates/GSB2.spec'
SB3 = 'galtemplates/GSB3.spec'
ES0 = 'galtemplates/Geso.spec'
Sbc = 'galtemplates/Gsbc.spec'
Scd = 'galtemplates/Gscd.spec'
lowz_galaxies = {'SB2':SB2, 'SB3':SB3, 'ES0':ES0, 'Sbc':Sbc, 'Scd':Scd}

def GalaxyTemplateColors(template, magfile, z0=4.0, z1=9.0, dz=0.1, ebmv=0.0, extlaw='xgal', filters=default_filters):
   # Calculate the colors as a function of redshift for all the filters 
   # All SEDs are normalized to have M_1500 = -21 mag (AB) before being 
   # redshifted. The magnitudes are output to magfile (it will be overwritten
   # if it already exists!)
   factory = FileSED.ZColorFactory(template, filters)
   # calculate all the magnitudes
   factory.colors_with_z(z0, z1, dz, ebmv=ebmv, extlaw=extlaw)
   f = open(magfile, 'wb')
   f.write('## template SED = %s\n' % template)
   f.write('## E(B-V) = %.2f\n' % ebmv)
   f.write('# 1 Z\n')
   for i in range(len(filters)):
      f.write('# %d %s\n' % (i+2, filters[i].name.split(',')[-1]))
   # now write the magnitudes
   for i in range(factory.num_z+1):
      f.write('%.2f ' % factory.zarray[i])
      for j in range(len(filters)):
         fname = filters[j].name.split(',')[-1]
         f.write('%.4f ' % factory.mags[fname][i])
      f.write('\n')
   f.close()

def LowZGalaxyColors(magfile, gtype='SB2', z0=0., z1=3.0, dz=0.05, ebmv=0.0, extlaw='xgal'):
   sedfile = lowz_galaxies[gtype]
   GalaxyTemplateColors(sedfile, magfile, z0=z0, z1=z1, dz=dz, ebmv=ebmv,
                        extlaw=extlaw)

class GalTemplateColors(object):
   # Given a galaxy template file and range of redshift, calculate the expected
   # observed magnitudes in the specified filters at each redshift (and their
   # colors).
   def __init__(self, template, refband=VJohnson, refmag=-21.0, ebmv=0., extlaw='xgal', filters=default_filters, cosmo=cosmo_def):
      self.template = template
      self.sp = S.FileSpectrum(template)
      if ebmv > 0:
         self.sp = self.sp * S.Extinction(ebmv, extlaw)
      if type(refband) == type(1500.):
         self.refband = S.Box(refband, 100.)
      else:
         self.refband = refband
      self.refmag = refmag
      self.filternames = [f.name for f in filters]
      self.kcorrs = {}
      for f in filters:
         self.kcorrs[f.name] = kcorr.KCorrect(self.sp, self.refband, f)
      if type(cosmo) == type([]):
         self.cosmo = cosmoclass.cosmoclass(H0=cosmo[0], omega_m=cosmo[1], omega_l=cosmo[2])
      else:
         self.cosmo = cosmo

   def calcMagnitudes(self, z):
      return np.array([self.refmag + self.cosmo.distmod(z) + self.kcorrs[fname](z) for fname in self.filternames])

   def writeMagnitudes(self, zarray, magfile=""):
      if not len(magfile):
         magfile = os.path.splitext(template)[0] + '.mag'
      f = open(magfile, 'wb')
      f.write('## %s\n' % self.template)
      f.write('# 1 Z\n')
      for i in range(len(self.filternames)):
         f.write('# %d %s\n' % (i+2, self.filternames[i].upper()))
      for j in range(len(zarray)):
         mags = self.calcMagnitudes(zarray[j])
         f.write('%f  ' % zarray[j])
         f.write('  '.join([str(m) for m in mags]) + '\n')
      f.close()


class PlotGalaxyColors(object):
   # Plot galaxy colors for a given SED whose magnitudes (as a function of z)
   # in different bands are given in magfile.
   def __init__(self, magfile):
      c = sextractor(magfile)
      for col in c._colnames:
         setattr(self, col, getattr(c, col))

   def get_color(self, b1, b2):
      assert hasattr(self, b1.lower()), "Magnitudes in %s does not exist." % b1.upper()
      assert hasattr(self, b2.lower()), "Magnitudes in %s does not exist." % b2.upper()
      mag1 = getattr(self, b1.lower())
      mag2 = getattr(self, b2.lower())
      return mag1 - mag2

   def plot_colors(self, color1, color2, z0, z1, ax=None, **plt_kwargs):
      # Plot color tracks from z0 to z1 (redshift step is determined by 
      # whatever self.magfile has).
      # color1 (and color2) is a list of names (e.g., ['f435w','f606w']) of filters in the 
      # format of [band1, band2], so the color is band1 - band2.
      # color1 appears at y and color2 appears at x.
      # All the other keyword arguments are passed to plt.plot
      assert z1 > z0, "z1 has to be greater than z0."
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      c1 = self.get_color(*color1)
      c2 = self.get_color(*color2)
      z0 = np.maximum(z0, self.z.min())
      z1 = np.minimum(z1, self.z.max())
      zrange = ((self.z >= z0) & (self.z < z1))
      ax.plot(c2[zrange], c1[zrange], **plt_kwargs)
      ax.set_xlabel('%s - %s' % (color2[0].upper(), color2[1].upper()))
      ax.set_ylabel('%s - %s' % (color1[0].upper(), color1[1].upper()))
      return ax

   def mark_redshifts(self, zmarks, color1, color2, ax=None, decimal=1, tol=1.e-4, **scatter_kwargs):
      # zmarks are redshifts to be marked; has to be in ascending order, and
      # preferable in regular steps
      zindices = []
      zmarked = []
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      for i in range(len(zmarks)):
         # find the index in self.z that matches zmarks[i]
         j = np.arange(len(self.z))[np.abs(self.z - zmarks[i])<=tol]
         if len(j) == 1:
            zindices += [j[0]]
            zmarked += [zmarks[i]]
      zindices = np.array(zindices)
      zmarked = np.array(zmarked)
      c1 = self.get_color(*color1)
      c2 = self.get_color(*color2)
      c1_marked = c1.take(zindices)
      c2_marked = c2.take(zindices)
      zstr = "%%.%df" % decimal
      z0 = zstr % zmarks[0]
      z1 = zstr % zmarks[-1]
      label = "z=%s ... %s" % (z0, z1)
      ax.scatter(c2_marked, c1_marked, label=label, **scatter_kwargs)
      # annotate the markers?
      return ax

class PlotLowZColors(PlotGalaxyColors):
   def __init__(self, magfile):
      super(PlotLowZColors, self).__init__(magfile)

   def plot_colors(self, color1, color2, z0=0., z1=3., ax=None, **plt_kwargs):
      ax=super(PlotLowZColors, self).plot_colors(color1, color2, z0=z0, z1=z1,
                                              ax=ax, **plt_kwargs)
      return ax
