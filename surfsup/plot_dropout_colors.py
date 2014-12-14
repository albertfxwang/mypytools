#!/usr/bin/env python

import numpy as np
from pygoods import Ftable, sextractor
import os, sys, glob
# from ImageTools import cutouts
from sedtools import galtemp_colors
from sedtools.bands import filters
import pysynphot as S
import matplotlib.pyplot as plt
import starlib
from bivariate2.plotting import plot_colors
from dropout_selection import lbg_colorcrit as lcc
# from pmscolors import pmscolors
import itertools

# f555w = S.ObsBandpass('acs,wfc2,f555w')
# f814w = S.ObsBandpass('acs,wfc2,f814w')
# f850lp = S.ObsBandpass('acs,wfc2,f850lp')
# f098m = S.ObsBandpass('wfc3,ir,f098m')
# f110w = S.ObsBandpass('wfc3,ir,f110w')
# f125w = S.ObsBandpass('wfc3,ir,f125w')
# f160w = S.ObsBandpass('wfc3,ir,f160w')
# filters = [f555w, f814w, f110w, f160w]
# filters_dic = {'f555w':f555w, 'f814w':f814w, 'f110w': f110w, 'f160w':f160w}
# zdrops_temp = 'lbgtemp_z6_0.4Zsolar_constSFH_100Myr.sed'
# where most of the galaxy template SEDs live
galtemp_dir = "/Users/khuang/Dropbox/codes/mypytools/sedtools/galtemplates"
LBG_template = os.path.join(galtemp_dir, 
                            "lbgtemp_z6_0.4Zsolar_constSFH_100Myr.sed")
magfiles_default = ['zdrops_macs0454_colors.dat', 
                  'zdrops_macs0454_colors_ext0.1.dat',
                  'zdrops_macs0454_colors_ext0.2.dat']
labels_default = ['E(B-V)=0', 'E(B-V)=0.1', 'E(B-V)=0.2']
sedtools = '/Users/khuang/Dropbox/codes/mypytools/sedtools'
lowz_magfiles = {'SB2': 'SB2colors_ext0.0_HST.dat',
               'SB3': 'SB3colors_ext0.0_HST.dat',
               'ES0': 'ES0colors_ext0.0_HST.dat',
               'Sbc': 'Sbccolors_ext0.0_HST.dat',
               'Scd': 'Scdcolors_ext0.0_HST.dat'}
starcolors = 'star_colors.dat'
markers_def = ['o', 'v', '^', '>', '<', '8', 's', 'p', 'D', 'h', 'd']
# p = pmscolors()
colors_def = ['Crimson','Green','DodgerBlue','HotPink','Lime','LightSkyBlue',
               'OrangeRed','SeaGreen','SteelBlue','Tomato','Aquamarine',
               'DarkOrchid']

# Try using plot_colors.LBGPlotColorFactory
class LBGPlotFactory(plot_colors.LBGColorPlotFactory):
   """
   A class that makes dropout-selection color-color plots.
   """
   def __init__(self, color1, color2, zcenter, magfiles=[], labels=[], galtemplate=LBG_template, ebmv=[0.,0.1,0.2,0.3], drop='f814w', droplabel='MACS0454', starcolors=plot_colors.starcolors_default):
      """
      color1, color2 : The filters that constitute each color, like 
                       color1 = ['f814w', 'f105w']
      zcenter        : The nominal central redshift for the target LBG sample.
      drop, droplabel: Two keywords for searching a specific set of color cuts
                       in LBG_criteria_data.yml
      starcolors     : A text file containing the expected colors of stars 
                       from the Pickles library; it may or may not have all 
                       the filters we want...
                       If starcolors == "", we will re-compute the colors 
                       within our custom filter set.
      magfiles       : Files containing the expected colors of galaxy
                       templates as a function of redshift. If magfiles == [],
                       we will re-compute the colors within our custom filter 
                       set.
      labels         : The label for each galaxy SED in magfiles. The length
                       of magfiles should equal the length of labels.
      """
      # magfiles = ['zdrops_macs0454_colors.dat',
      #             'zdrops_macs0454_colors_ext0.1.dat',
      #             'zdrops_macs0454_colors_ext0.2.dat']
      # labels = ['E(B-V)=0','E(B-V)=0.1','E(B-V)=0.2']
      assert len(magfiles) == len(labels)
      self.zcenter = zcenter
      self.color1 = map(lambda x: x.lower(), color1)
      self.color2 = map(lambda x: x.lower(), color2)
      self.drop = drop
      if not len(magfiles):
         for i in range(len(ebmv)):
            mf_name = '%s_drops_%s_colors_ext%.2f.dat' % (drop, droplabel, ebmv[i])
            self.calc_LBG_colors(galtemplate, mf_name, ebmv=ebmv[i])
            magfiles += [mf_name]
            labels += ['E(B-V)=%.2f' % ebmv[i]]
      self.zmarks = np.linspace(zcenter-1., zcenter+1., 5)
      print self.zmarks
      if not starcolors:
         raise NotImplementedError, "Please calculate the expected stellar colors first..."
      super(LBGPlotFactory, self).__init__(magfiles, labels, self.zcenter, 
                                          zmarks=self.zmarks, 
                                          starcolors=starcolors)
      # can also plot the dropout criteria
      self.lcc = lcc.colorcrit()
      self.lcc = self.lcc(drop, droplabel)      

   def calc_LBG_colors(self, galtemplate, magfile, zrange=None, dz=0.05, ebmv=0.):
      if zrange == None:
         z0 = self.zcenter - 2.
         z1 = self.zcenter + 2.
      else:
         z0, z1 = zrange
      lbg_filters = [filters[self.color1[0]], filters[self.color1[1]]]
      if self.color2[0] != self.color1[1]:
         lbg_filters += [filters[self.color2[0]]]
      lbg_filters += [filters[self.color2[1]]]
      galtemp_colors.GalaxyTemplateColors(galtemplate, magfile, 
                     z0=z0, z1=z1, dz=dz, ebmv=ebmv, filters=lbg_filters)

   def plot_tracks(self, ax=None, zannotate=-1, sizeannotate='xx-large', lowz_label=False, zlabel=True, offset=[-0.5,0.5], **plt_kwargs):
      if zannotate < 0:
         zannotate = self.zcenter
      ax = super(LBGPlotFactory, self).plot_tracks(self.color1, self.color2, 
                  ax=ax, zmarks=self.zmarks, zannotate=zannotate, 
                  sizeannotate=sizeannotate, zlabel=zlabel, offset=offset, 
                  **plt_kwargs)
      ax = LBGPlotFactory.plot_lowz(self, self.color1, self.color2, ax=ax,
                                    label=lowz_label, **plt_kwargs)
      return ax

   def plot_all(self, field_name, figsize=(10,8), critlabel="", critlabelsize='xx-large', xmin=-1, xmax=3, ymin=-1., ymax=9., fc_cc='blue', alpha_cc=0.2, sizeannotate='xx-large', zannotate=None, lowz_label=False, zlabel=True, offset=[-0.5,0.5], **plt_kwargs):
      fig = plt.figure(figsize=figsize)
      ax = fig.add_subplot(111)
      ax.set_xlim(xmin, xmax)
      ax.set_ylim(ymin, ymax)
      ax = self.plot_colorcrit(self.lcc, ax, fc=fc_cc, alpha=alpha_cc)
      ax = self.plot_tracks(zannotate=zannotate, sizeannotate=sizeannotate, 
                            ax=ax, lowz_label=lowz_label, zlabel=zlabel, 
                            offset=offset, **plt_kwargs)
      ax.set_xlim(xmin, xmax)
      ax.set_ylim(ymin, ymax)
      ax.set_title('%s-dropouts @ z~%.1f in %s' % (self.drop.upper(), self.zcenter, field_name))
      ax.set_xlabel(ax.get_xlabel(), size='xx-large')
      ax.set_ylabel(ax.get_ylabel(), size='xx-large')
      if len(critlabel):
         ax.text(0.05, 0.95, critlabel, size=critlabelsize, 
                 ha='left', va='top', 
                 transform=ax.transAxes, 
                 bbox=dict(boxstyle='round', fc='none'))
      return ax

   def plot_galaxy(self, ax, c, objid, magform='iso', **ebar_kwargs):
      match = (c.number==objid)
      mags = []
      magerrs = []
      bands = self.color1 + self.color2
      for i in range(len(bands)):
         b = bands[i]
         mag = getattr(c,'%s_mag_%s_1sig'%(b,magform))[match][0]
         magerr = getattr(c,'%s_magerr_%s'%(b,magform))[match][0]
         flux = getattr(c,'%s_flux_%s'%(b,magform))[match][0]
         fluxerr = getattr(c,'%s_fluxerr_%s_scaled'%(b,magform))[match][0]
         if (flux / fluxerr) <= 1.:
            magerr = -1.0
         mags += [mag]
         magerrs += [magerr]
      color1 = mags[0] - mags[1]
      if (magerrs[0] < 0) or np.abs(mags[0]-magerrs[0])<=1.e-4:
         color1_err = -1.0
      else:
         color1_err = np.sqrt(magerrs[0]**2 + magerrs[1]**2)
      color2 = mags[2] - mags[3]
      if (magerrs[2] < 0) or np.abs(mags[2]-magerrs[2])<=1.e-4:
         color2_err = -1.0
      else:
         color2_err = np.sqrt(magerrs[2]**2 + magerrs[3]**2)
      print "Object ID # = %d" % objid
      print "%s-%s = %.2f +/- %.2f" % (bands[0], bands[1], color1, color1_err)
      print "%s-%s = %.2f +/- %.2f" % (bands[2], bands[3], color2, color2_err)
      ax=super(LBGPlotFactory, self).plot_galaxy(ax, color1, color1_err,
                                                color2, color2_err, 
                                                **ebar_kwargs)
      return ax

   def plot_all_galaxies(self, ax, c, magform='iso', color='0.5', size=4, marker='o'):
      """
      Plot all galaxies in the field as tiny gray dots (or some other marker).
      Input c should be a hst_catalog.HSTCatalog instance.
      """
      color1_all = [c.calc_color(x, self.color1[0], self.color1[1]) for x in c.number]
      color1_all = np.array(color1_all)[:,0]
      color2_all = [c.calc_color(x, self.color2[0], self.color2[1]) for x in c.number]
      color2_all = np.array(color2_all)[:,0]
      ax.scatter(color2_all, color1_all, color=color, s=size, marker=marker)
      return ax

   def plot_candidates(self, ax, c, objids, ms_show=14, objnames=[], magform='iso', markers=[], ncol_legend=1,objnames2show=[], marker_noshow='s', ms_noshow=8, color_noshow='black', **ebar_kwargs):
      # Re-start iterating markers... so markers are consistent for the same 
      # objects
      if len(markers):
        assert len(markers) == len(objids)
        markers_iter = itertools.cycle(markers)
      else:
        markers_iter = itertools.cycle(markers_def)
      colors_iter = itertools.cycle(colors_def)
      if not len(objnames):
         objnames = map(lambda x: str(x), objids)
      for i in range(len(objnames)):
         objLabel = ""
         marker = marker_noshow
         color = color_noshow
         ms = ms_noshow
         mew = 1
         ebar_kwargs['elinewidth'] = 1
         if objnames[i] in objnames2show:
            objLabel = objnames[i]
            marker = markers_iter.next()
            color = colors_iter.next()
            ms = ms_show
            mew = 2
            ebar_kwargs['elinewidth'] = 2
         ax = self.plot_galaxy(ax, c, objids[i], magform=magform, 
                               label=objLabel, ms=ms, mew=mew,
                               color=color, 
                               fmt=marker, **ebar_kwargs)
      ax.legend(loc=4, fontsize='medium', ncol=ncol_legend)
      plt.draw()
      return ax   

def stars_colors(magfile=starcolors, normband='f160w'):
  # calculate colors of stars
  starlib.calcmag(magfile, filters=filters_dic, normband=normband)


# def plot_zdrops_colors(magfiles=magfiles_default, labels=labels_default, color1=['f814w','f110w'], color2=['f110w','f160w'], galtemplate=zdrops_temp, z0=4.0, z1=9.0, zmarks=[6.0, 6.5, 7.0, 7.5], legendloc=4, ax=None):
#    # Plot color tracks
#    # color1 (and same with color2) should be like ['f435w', 'f814w'] (two 
#    # filter names in strings); this means the color is f435w - f814w
#    if ax==None:
#       fig = plt.figure()
#       ax = fig.add_subplot(111)
#    for i in range(len(magfiles)):
#       P = galtemp_colors.PlotGalaxyColors(magfiles[i])
#       ax = P.plot_colors(color1, color2, z0=z0, z1=z1, label=labels[i], ax=ax)
#       if i == 0:
#          # also mark certain redshift steps
#          ax = P.mark_redshifts(zmarks, color1, color2, ax=ax, 
#                                facecolor='none', edgecolor='red', s=10**2, 
#                                marker='o')
#    # Now plot low-z galaxies
#    for g in lowz_magfiles.keys():
#       P = galtemp_colors.PlotLowZColors(sedtools+'/'+lowz_magfiles[g])
#       ax = P.plot_colors(color1, color2, label=g, ax=ax, ls='--')
#    # Now plot stars
#    pstars = starlib.starlib(starcolors)
#    c1 = pstars.colors(*color1)
#    c2 = pstars.colors(*color2)
#    ax.scatter(c2, c1, marker='*', s=8**2, label='stars', facecolor='yellow')

#    ax.legend(loc=legendloc, scatterpoints=1)
#    ax.set_ylim(-1.,6.)
#    ax.set_xlim(xmax=8.)
#    plt.title('Colors for z-dropouts at z~7', size=22)


