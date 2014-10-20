#!/usr/bin/env python

from numpy import *
from pygoods import *
import pysynphot as S
import magredshift as MR
from pysynphot import Extinction
import myutil
import arrayspec as A
import simutil

Lsun = 3.826e33  # erg / sec
pc_cm = 30.857e17  # 1 parsec in cm
#uni1500 = S.FileBandpass('uni1500.bp')
uni1500 = S.Box(1500.,100.)
iband = S.ObsBandpass('acs,wfc2,f775w')
zband = S.ObsBandpass('acs,wfc2,f850lp')
#sp10 = S.FileSpectrum('mylbg_sfr10.sed')


def simkcorr(sp, absmag_rest, restband, obsband, z, lya_ew):
   # defined as m_1 = M_2 + kcorr
   # m_1: apparent magnitude in observed band
   # M_2: absolute magnitude in rest-frame band
   # assumes BC03 model: model flux unit is in L_solar/angstrom
   #sp = S.FileSpectrum(spec)
   tenpc = 4.*pi*(10.*pc_cm)**2  # 10 parsec in centimeter
   absmag = A.ABmag(sp,restband)
   spn = sp * 10.**((absmag_rest-absmag)/-2.5)  # normalize sp to absmag_rest
   spn = spn * tenpc  # get luminosity
   # Now redshift spn to z
   obsmag,spec_at,spec_r = MR.mag_redshift(spn,z,obsband,"abmag")
   lyalam = 1216. * (1. + z)  # observed wavelength of Lya line
   spec_lya = simutil.add_lyman_alpha(spec_at, lyalam, lya_ew)
   obsmag = A.ABmag(spec_lya, obsband)
   kc = obsmag - absmag_rest
   return obsmag, kc


def mconvertfile(sp, zlo, zhi, dz, mcfile, restband=uni1500, obsband=iband):
   zarr = arange(zlo, zhi+dz, dz); print len(zarr)
   dm = zeros(len(zarr))
   for i in range(len(zarr)):
      obsmag, kc = simkcorr(sp, -21.0, restband, obsband, zarr[i])
      dm[i] = kc
   f = open(mcfile,'w')
   f.write('# %10s %15s\n' % ('z', 'dm'))
   for i in range(len(dm)):
      f.write('  %10f %15f\n' % (zarr[i], dm[i]))
   f.close()
   return 0 


def make_mconvert_1500_i(zlo, zhi, dz, mcfile, ebmv=0.15):
   calz = S.Extinction(ebmv, 'xgal')  # dust obscuration
   sp_ext = sp10 * calz  # dust obscured SED
   mconvertfile(sp_ext, zlo, zhi, dz, mcfile, restband=uni1500, obsband=iband)
   return 0


def make_mconvert_1500_z(zlo, zhi, dz, mcfile, ebmv=0.15):
   calz = S.Extinction(ebmv, 'xgal')  # dust obscuration
   sp_ext = sp10 * calz  # dust obscured SED
   mconvertfile(sp_ext, zlo, zhi, dz, mcfile, restband=uni1500, obsband=zband)
   return 0 


