#!/usr/bin/env python

from numpy import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from pygoods import Ftable, sextractor
from dropout_selection import lbg_colorcrit as lcc
import hst_zeropoints
from my_mplfonts import Helvetica

"""
Define parent class plot_dropouts.
"""

class plot_dropouts(object):
   def __init__(self,catalog):
      # matplotlib keyword arguments
      self.c = Ftable(catalog)
      # Dropout catalog as a FITS table
      self.histkw = {'histtype':'step','lw':1.5}
      self.errorbarkw = {'ms':6,'capsize':6,'mec':'blue','ecolor':'blue'}

   def plot_photoz_gfmag(self, field='udf', gfband='f606w', 
                         dropname='U-dropouts'):
      field = field.lower()
      field_crit = getattr(self.c, '%s_fit' % field.lower())
      magout = getattr(self.c, '%s_magout_gf'%gfband)[field_crit==True]
      photo_z = self.c.photo_z[field_crit==True]
      photo_z_upper_68 = self.c.photo_z_upper_68[field_crit==True]
      photo_z_plus = photo_z_upper_68 - photo_z
      photo_z_lower_68 = self.c.photo_z_lower_68[field_crit==True]
      photo_z_minus = photo_z - photo_z_lower_68
      #print "photo_z_minus", photo_z_minus
      fig = plt.figure()
      ax = fig.add_subplot(1,1,1)
      ax.errorbar(magout, photo_z, yerr=[photo_z_minus, photo_z_plus],
                  fmt='s', capsize=5)
      spec_z_crit = (field_crit==True) & (in1d(self.c.spec_z_dq, (1,2)))
      if sum(spec_z_crit) > 0:
         spec_z = self.c.spec_z[spec_z_crit==True]
         magout_specz = getattr(self.c, '%s_magout_gf'%gfband)[spec_z_crit==True]
         ax.scatter(magout_specz, spec_z, marker='*', s=10**2, 
                    color='orange')
      ax.set_xlabel('%s GALFIT mag'%gfband.upper(), size=18)
      ax.set_ylabel('Photo-z (Bayesian)', size=18)
      ax.set_title('%s in %s' % (dropname,field.upper()), size=22)

class plot_goods_acs_v2(object):
   """
   Plot the GOODS/ACS v2 catalogs.
   """
   def __init__(self, cb, cv, ci, cz, field='n', drop='b', key='g04'):
      """
      Read the catalogs in all four filters.
      """
      self.c = {}
      self.bands = ['f435w','f606w','f775w','f850lp']
      self.c['f435w'] = cb
      self.c['f606w'] = cv
      self.c['f775w'] = ci
      self.c['f850lp'] = cz
      if field in ['n', 's']:
         self.zeropoints = hst_zeropoints.zpt_goodsv2
      elif field == 'udf':
         self.zeropoints = hst_zeropoints.zpt_udf
      self.magiso_1sig = {}
      self.magauto_1sig = {}
      self.drop = drop
      for b in self.bands:
         # calculate 1-sigma FLUX_ISO
         flux_iso = self.c[b].flux_iso
         fluxerr_iso = self.c[b].fluxerr_iso
         flux_iso_1sig = where(flux_iso/fluxerr_iso > 1., flux_iso, 
                               fluxerr_iso)
         mag_iso_1sig = self.zeropoints[b] - 2.5 * log10(flux_iso_1sig)
         self.magiso_1sig[b] = mag_iso_1sig
         # calculate 1-sigma FLUX_AUTO
         flux_auto = self.c[b].flux_auto
         fluxerr_auto = self.c[b].fluxerr_auto
         flux_auto_1sig = where(flux_auto/fluxerr_auto > 1., flux_auto,
                                fluxerr_auto)
         mag_auto_1sig = self.zeropoints[b] - 2.5 * log10(flux_auto_1sig)
         self.magauto_1sig[b] = mag_auto_1sig
      self.lcc_dropout = lcc.colorcrit()
      self.lcc_dropout(drop, key)

   def plot_2colors(self, color1, color2, colorname1, colorname2, ax=None,
                    xlims=[-2.,4.], ylims=[-2.,5.], vmin=20., vmax=28.2):
      cb = self.c['f435w']
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      self.lcc_dropout.select(color1, color2)
      # robust colors but not dropouts
      self.detect_nd = (self.lcc_dropout.crit==False) & (self.detect==True) & \
         (self.uplim==False)
      # robust colors and dropouts
      self.detect_d = (self.lcc_dropout.crit==True) & (self.detect==True) & \
         (self.uplim==False)
      # 1-sig colors but not dropouts
      self.uplim_nd = (self.lcc_dropout.crit==False) & (self.detect==True) & \
         (self.uplim==True)
      # 1-sig colors and dropouts
      self.uplim_d = (self.lcc_dropout.crit==True) & (self.detect==True) & \
         (self.uplim==True)
      color1_detect_nd = color1[self.detect_nd==True]
      color1_detect_d = color1[self.detect_d==True]
      color1_uplim_nd = color1[self.uplim_nd==True]
      color1_uplim_d = color1[self.uplim_d==True]
      color2_detect_nd = color2[self.detect_nd==True]
      color2_detect_d = color2[self.detect_d==True]
      color2_uplim_nd = color2[self.uplim_nd==True]
      color2_uplim_d = color2[self.uplim_d==True]
      zmag = self.magauto_1sig['f850lp']
      zmag_detect_d = zmag[self.detect_d==True]
      zmag_uplim_d = zmag[self.uplim_d==True]

      # Plot the non-dropouts as gray
      cax1 = ax.scatter(color2_detect_nd, color1_detect_nd, marker='o', s=2**2, 
                 edgecolors='none', c='0.5')
      cax2 = ax.scatter(color2_uplim_nd, color1_uplim_nd, marker='^', s=2**2,
                 edgecolors='none', c='0.5')
      # Plot the dropouts color-coded by F850LP magnitude
      cax3 = ax.scatter(color2_detect_d, color1_detect_d, marker='o', s=5**2, 
                 edgecolors='none', c=zmag_detect_d, vmin=vmin, vmax=vmax)
      cax4 = ax.scatter(color2_uplim_d, color1_uplim_d, marker='^', s=5**2,
                 edgecolors='none', c=zmag_uplim_d, vmin=vmin, vmax=vmax)
      cbar = plt.colorbar(cax3)
      cbar.set_label('F850LP magnitude', font_properties=Helvetica(13))
      cbar.set_ticks(arange(vmin, vmax+0.2))
      self.lcc_dropout.plotcrit(xlims[0], ylims[1], ax=ax, ls='--', lw=2.0)
      ax.set_xlabel(colorname2, font_properties=Helvetica(16))
      ax.set_ylabel(colorname1, font_properties=Helvetica(16))
      xticks = arange(xlims[0], xlims[1]+1.)
      yticks = arange(ylims[0], ylims[1]+1.)
      ax.set_xticks(xticks)
      ax.set_yticks(yticks)
      ax.set_xticklabels(xticks, map(lambda x:'%d'%x, xticks),
                         font_properties=Helvetica(13))
      ax.set_yticklabels(yticks, map(lambda y:'%d'%y, yticks),
                         font_properties=Helvetica(13))
      ax.set_xlim(*xlims)
      ax.set_ylim(*ylims)
      ax.set_title('%s-dropouts' % self.drop.upper(), 
                   font_properties=Helvetica(20))

class plot_bdrops_goods_acs_v2(plot_goods_acs_v2):
   def __init__(self, cb, cv, ci, cz, field='n', key='g04'):
      plot_goods_acs_v2.__init__(self, cb, cv, ci, cz, field=field, drop='b',
                                 key=key)
      self.detect = (cv.flux_auto>0) & (ci.flux_auto>0) & (cz.flux_auto>0)
      self.detect = (self.detect & (cz.flux_auto/cz.fluxerr_auto >= 5.0))
      self.detect = (self.detect & (cz.mag_auto <= 30.0))
      self.detect = (self.detect & (cv.flux_auto/cv.fluxerr_auto >= 2.0))
      self.uplim = (cb.flux_iso/cb.fluxerr_iso <= 1.0)
      print "sum(self.detect)", sum(self.detect)
      print "sum(self.uplim)", sum(self.uplim)

   def plot_bmv_vmz(self, xlims=[-1.,4.], ylims=[-1.,5.], vmin=20., vmax=28.2,
                    ax=None):
      bmv = (self.magiso_1sig['f435w']-self.magiso_1sig['f606w'])
      vmz = (self.magiso_1sig['f606w']-self.magiso_1sig['f850lp'])
      self.plot_2colors(bmv, vmz, 'B - V', 'V - z', xlims=xlims, 
                        ylims=ylims, vmin=vmin, vmax=vmax, ax=ax)

class plot_vdrops_goods_acs_v2(plot_goods_acs_v2):
   def __init__(self, cb, cv, ci, cz, field='n', key='g04'):
      plot_goods_acs_v2.__init__(self, cb, cv, ci, cz, field=field, drop='v',
                                 key=key)
      self.detect = (ci.flux_auto>0) & (cz.flux_auto>0)
      self.detect = (self.detect & (cz.flux_auto/cz.fluxerr_auto >= 5.0))
      self.detect = (self.detect & (cz.mag_auto <= 30.0))
      self.detect = (self.detect & (ci.flux_auto/ci.fluxerr_auto >= 2.0))
      self.uplim = (cv.flux_iso/cv.fluxerr_iso <= 1.0)
      print "sum(self.detect)", sum(self.detect)
      print "sum(self.uplim)", sum(self.uplim)

   def plot_vmi_imz(self, ax=None, xlims=[-3.,4.], ylims=[-1.,6.], vmin=20.,
                    vmax=28.2):
      vmi = (self.magiso_1sig['f606w']-self.magiso_1sig['f775w'])
      imz = (self.magiso_1sig['f775w']-self.magiso_1sig['f850lp'])
      self.plot_2colors(vmi, imz, 'V - i', 'i - z', xlims=xlims,
                        ylims=ylims, vmin=vmin, vmax=vmax, ax=ax)
