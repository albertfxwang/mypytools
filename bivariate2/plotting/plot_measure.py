#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
from pygoods import Ftable, sextractor
from my_mplfonts import Helvetica
from bivariate2 import galaxy_samples as gs
from plotting import plot_catalog
import os, sys

"""
Plot GALFIT/SExtractor measurements of galaxy samples.
"""
# Set some default plot properties
mag_lim = [22.0, 29.0]
re_lim = [0.01, 3.]
# Set default marker properties
galfit_fc = 'blue'  # facecolor for GALFIT points
galfit_ec = 'black'  # edgecolor for GALFIT points
galfit_ms = 6**2  # marker size for GALFIT points
sex_fc = '0.5'    # facecolor for SExtractor points
sex_ec = '0.5'  # edgecolor for SExtractor_points
sex_ms = 4**2     # marker size for SExtractor points
## U-dropouts
udrops_field_names = ['udf', 'deep+ers', 'wide']
udrops_galfit_fc = ['red', 'green', 'blue']
udrops_catalog_dir = '/Users/khuang/Dropbox/Research/bivariate/udrops_sample'
udrops_catalogs = ['gds_udrops_hudf_30mas_140313_galfit.fits',
                  'gds_udrops_deep+ers_30mas_140313_galfit.fits',
                  'gds_udrops_wide_30mas_140313_galfit.fits']
udrops_gfmag_lims = [28.5, 26.5, 26.5]
udrops_markers = ['o', '^', 's']
## B-dropouts
bdrops_field_names = ['udf', 'goods']
bdrops_galfit_fc = ['red', 'blue']
bdrops_catalog_dir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/bdrops_fitting/catalogs'
bdrops_catalogs = ['bdrops_udf_gf_v3.fits',
                   'bdrops_gf_v3.fits']
bdrops_gfmag_lims = [28.5, 26.5]
bdrops_markers = ['o', 's']
## V-dropouts
vdrops_field_names = ['udf', 'goods']
vdrops_galfit_fc = ['red', 'blue']
vdrops_catalog_dir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/vdrops_fitting/catalogs'
vdrops_catalogs = ['vdrops_udf_gf_v3.fits',
                   'vdrops_gf_v3.fits']
vdrops_gfmag_lims = [28.5, 26.5]
vdrops_markers = ['o', 's']         
## i-dropouts
idrops_field_names = ['udf', 'deep+ers', 'wide']
idrops_galfit_fc = ['red', 'green', 'blue']
idrops_catalog_dir = '/Users/khuang/Dropbox/Research/bivariate/idrops_sample'
idrops_catalogs = ['idrops_hudf_gf_f125w_140228.fits',
                      'idrops_gds_deep+ers_gf_f125w_140228.fits',
                      'idrops_gds_wide_gf_f125w_140228.fits']
idrops_gfmag_lims = [28.0, 27.0, 26.5]  
idrops_markers = ['o', '^', 's']    

class PlotOneField(plot_catalog.PlotCatalog):
   def __init__(self, catalog, field_name, filter_name, format='FITS', good_vflags=[0,1,3,12]):
      # Plot columns for a single field
      # Accepts multiple catalogs representing different fields at once.
      assert (format in ['FITS','SEX']), "Only FITS tables and SExtractor catalogs are supported."
      if format == 'FITS':
         readcat = Ftable
      else:
         readcat = sextractor
      self.c = readcat(catalog)
      # self.sample_name = sample_name
      self.field_name = field_name.lower()
      self.filter_name = filter_name.lower()
      self.titlesize = 20
      self.axlabelsize = 22
      self.ticklabelsize = 16
      vflag = self.c.vflag
      self.vflag_sample = np.in1d(vflag, good_vflags)
      self.galfit_goodfits = self.goodfits() & np.in1d(vflag, good_vflags)
      self.galfit_badfits = (self.vflag_sample==True) & \
                            (self.galfit_goodfits==False)

   def goodfits(self, reerr_lim=0.6, chi2nu_lim=5.0, re_col='re_out_arcsec', reerr_col='re_out_err_arcsec', chi2nu_col='chi2nu'):
      # determine good GALFIT fits
      re_out = self.get_col('%s_%s' % (self.filter_name, re_col))
      re_err = self.get_col('%s_%s' % (self.filter_name, reerr_col))
      chi2nu = self.get_col('%s_%s' % (self.filter_name, chi2nu_col))
      magflag = self.get_col('%s_magflag' % self.filter_name)
      reflag = self.get_col('%s_reflag' % self.filter_name)
      nflag = self.get_col('%s_nflag' % self.filter_name)
      good = (re_err / re_out <= reerr_lim) & (re_out > 0.) & (re_err > 0.)
      good = good & (chi2nu <= chi2nu_lim) & (magflag == 0) & (reflag == 0)
      good = good & (nflag == 0)
      return good

   def get_galfit_col(self, col_base, condition=None):
      # get GALFIT-measured columns
      # By default, only return the good ones
      if condition == None:
         condition = self.galfit_goodfits
      col = self.get_col('%s_%s' % (self.filter_name, col_base),
                         condition=condition)
      return col

   def get_sex_col(self, col_base, condition=None):
      # get SExtractor-measured columns
      # By default, only return those with POOR GALFIT fits
      if condition == None:
         # condition = np.logical_not(self.galfit_goodfits)
         condition = self.galfit_badfits
      col = self.get_col('%s_%s' % (self.filter_name, col_base),
                         condition=condition)
      return col

   def scatter_galfit_mag_logre(self, mag_col='mag_out', re_col='re_out_arcsec', condition=None, ax=None, label='GALFIT', **scatterkw):
      if condition == None:
         condition = self.galfit_goodfits
      mag_colname = '%s_%s' % (self.filter_name, mag_col)
      re_colname = '%s_%s' % (self.filter_name, re_col)
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      ax.set_yscale('log') # set Re to log scale
      ax = self.scatter_xy_cols(mag_colname, re_colname, 
                           condition=condition, ax=ax, label=label, 
                           **scatterkw)
      plt.draw()
      return ax

   def scatter_sex_mag_logre(self, mag_col='mag_auto', re_col='flux_radius_2_arcsec', condition=None, ax=None, label='SExtractor', **scatterkw):
      if condition == None:
         condition = np.logical_not(self.galfit_goodfits)
      mag_colname = '%s_%s' % (self.filter_name, mag_col)
      re_colname = '%s_%s' % (self.filter_name, re_col)
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      ax.set_yscale('log') # set Re to log scale
      ax = self.scatter_xy_cols(mag_colname, re_colname, 
                           condition=condition, ax=ax, label=label, 
                           **scatterkw)
      plt.draw()
      return ax

   def set_labels(self, ax, xname, yname, xlim=None, ylim=None):
      # Set label names
      ax.set_xlabel(xname, size=self.axlabelsize)
      ax.set_ylabel(yname, size=self.axlabelsize)
      if xlim != None:
         ax.set_xlim(xlim)
      if ylim != None:
         ax.set_ylim(ylim)
      plt.draw()
      return ax

   def add_gfmag_lim(self, ax, mag_lim, **plot_kwargs):
      ## Add a vertical line for GALFIT magnitude limits
      ylims = ax.get_ylim()
      ax.plot([mag_lim, mag_lim], ylims, **plot_kwargs)
      return ax

class PlotFactory(object):
   def __init__(self, catalogs, field_names, filter_name, sample="", colors=None, catalog_dir='.', gfmag_lims=[], **pof_kwargs):
      assert len(catalogs) == len(field_names)
      self.plots = {}
      self.field_names = field_names
      self.filter_name = filter_name
      self.sample = sample
      self.colors = colors
      if len(gfmag_lims) == 0:
         gfmag_lims = [-1] * len(field_names)
      self.gfmag_lims = gfmag_lims
      for i in range(len(field_names)):
         f = field_names[i]
         cat = os.path.join(catalog_dir, catalogs[i])
         self.plots[f] = PlotOneField(cat, f, filter_name, **pof_kwargs)

   def plot_grid(self, figsize=(12,4), nrows=1, suptitle="", axes_pad=0.05):
      # Make a grid of plots, with each panel showing one field
      fig = plt.figure(figsize=figsize)
      ncols = len(self.field_names) / nrows
      if len(self.field_names) % nrows != 0:
         ncols += 1
      print "nrows, ncols = ", nrows, ncols
      grid = ImageGrid(fig, (0.1,0.15,0.8,0.8), nrows_ncols=(nrows, ncols), 
                       axes_pad=axes_pad, aspect=False)
      for i in range(len(self.field_names)):
         ax_i = grid[i]
         f = self.field_names[i]
         p = self.plots[f]
         # First plot GALFIT points
         p.scatter_galfit_mag_logre(ax=ax_i, label="", marker='s', 
                                    facecolor=galfit_fc, edgecolor=galfit_ec,
                                    s=galfit_ms)
         p.scatter_sex_mag_logre(ax=ax_i, marker='o', facecolor=sex_fc, 
                                 edgecolor=sex_ec, s=sex_ms, label="")
         ax_i.text(0.05, 0.05, f.upper(), transform=ax_i.transAxes, 
                   va='bottom', ha='left', size='x-large',
                   bbox=dict(boxstyle='round', fc='none'))
         ax_i.set_xticklabels(np.arange(mag_lim[0],mag_lim[-1]).astype('int'))
         ax_i.set_xlabel('%s magnitude' % self.filter_name.upper())
         ax_i.set_ylabel('%s galaxy size [arcsec]' % self.filter_name.upper())
         if self.gfmag_lims[i] > 0:
            x = self.gfmag_lims[i]
            ax_i.plot([x,x], re_lim, ls='--', c='0.2', lw=2.0)
         ax_i.set_xlim(mag_lim)
         ax_i.set_ylim(re_lim)
      plt.draw()
      if len(suptitle):
         fig.suptitle(suptitle)
      return grid

   def plot_single(self, markers, ax=None, textsize=20):
      # Plot the measurements from all fields (for the same galaxy sample) 
      # in the same panel
      # Use different marker shapes to distinguish between GALFIT measurments 
      # in different fields
      if ax == None:
         fig = plt.figure()
         ax = fig.add_subplot(111)
      assert len(markers) == len(self.field_names)
      ax.set_ylim(0.005, 20.0)
      for i in range(len(markers)):
         m = markers[i]
         f = self.field_names[i]
         p = self.plots[f]
         ax = p.scatter_galfit_mag_logre(ax=ax, label=f.upper(), marker=m, 
                                    s=galfit_ms, facecolor=self.colors[i],
                                    edgecolor=self.colors[i])
         ax = p.scatter_sex_mag_logre(ax=ax, label="", marker='x', 
                                 edgecolor=sex_ec, s=sex_ms)
         ax = p.add_gfmag_lim(ax, self.gfmag_lims[i], lw=2.0, ls='--', 
                              color=self.colors[i])
      ax.set_xlabel("magnitude")
      ax.set_ylabel("angular size [arcsec]")
      ax.text(0.05, 0.95, "%s\n(%s)" % (self.sample, self.filter_name.upper()), 
              transform=ax.transAxes, ha='left', 
              va='top', size=textsize, 
              bbox=dict(boxstyle='round', facecolor='none'))
      ax.legend(loc=3)
      # Add text for each GALFIT magnitude limit
      ylims = ax.get_ylim()
      delta_logy = np.log10(ylims[1]) - np.log10(ylims[0])
      ## First, UDF-depth
      ax.text(self.gfmag_lims[0]-0.01, 3.0,
              # 10.**(np.log10(ylims[0])+0.95*delta_logy), 
              r'UDF$\rightarrow$', 
              color=self.colors[0], size=textsize,
              ha='right', va='center')
      ## Then, the shallower depths
      fn = r'%s$\rightarrow$' % self.field_names[1].upper()
      if (len(self.field_names) > 2):
         if (self.gfmag_lims[1] == self.gfmag_lims[2]):
            fn = r'%s+%s$\rightarrow$' % (self.field_names[1].upper(), self.field_names[2].upper())         
      ax.text(self.gfmag_lims[1]-0.01, 2.0,
         # 10.**(np.log10(ylims[0])+0.9*delta_logy),
         fn, size=textsize, 
         ha='right', va='center', color=self.colors[1])
      # If there is a third different depth (e.g., for i-dropouts)
      if (len(self.field_names) > 2):
         if (self.gfmag_lims[1] != self.gfmag_lims[2]):
            fn3 = r'%s$\rightarrow$' % self.field_names[2].upper()
            ax.text(self.gfmag_lims[2]-0.01, 1.0,
                    # 10.**(np.log10(ylims[0])+0.85*delta_logy),
                    fn3, size=textsize, 
                    ha='right', va='center', color=self.colors[2])
      ax.set_xlim(21., 29.5)
      plt.draw()
      return ax

class UdropsPlotFactory(PlotFactory):
   def __init__(self):
      super(UdropsPlotFactory, self).__init__(udrops_catalogs, 
                              udrops_field_names,
                              'f606w',
                              sample='z~3.2',
                              colors=udrops_galfit_fc,
                              catalog_dir=udrops_catalog_dir,
                              gfmag_lims=udrops_gfmag_lims)

   def plot_single(self, ax=None, textsize=20):
      ax = super(UdropsPlotFactory, self).plot_single(udrops_markers, ax=ax, 
                                                 textsize=textsize)
      return ax

class BdropsPlotFactory(PlotFactory):
   def __init__(self):
      super(BdropsPlotFactory, self).__init__(bdrops_catalogs,
                              bdrops_field_names,
                              'f775w',
                              sample='z~4',
                              colors=bdrops_galfit_fc,
                              catalog_dir=bdrops_catalog_dir,
                              gfmag_lims=bdrops_gfmag_lims)

   def plot_single(self, ax=None, textsize=20):
      ax = super(BdropsPlotFactory, self).plot_single(bdrops_markers, ax=ax,
                                                 textsize=textsize)
      return ax

class VdropsPlotFactory(PlotFactory):
   def __init__(self):
      super(VdropsPlotFactory, self).__init__(vdrops_catalogs,
                              vdrops_field_names,
                              'f850lp',
                              sample='z~5',
                              colors=vdrops_galfit_fc,
                              catalog_dir=vdrops_catalog_dir,
                              gfmag_lims=vdrops_gfmag_lims)

   def plot_single(self, ax=None, textsize=20):
      ax = super(VdropsPlotFactory, self).plot_single(vdrops_markers, ax=ax,
                                                 textsize=textsize)
      return ax

class idropsPlotFactory(PlotFactory):
   def __init__(self):
      super(idropsPlotFactory, self).__init__(idrops_catalogs,
                              idrops_field_names,
                              'f125w',
                              sample='z~6',
                              colors=idrops_galfit_fc,
                              catalog_dir=idrops_catalog_dir,
                              gfmag_lims=idrops_gfmag_lims)

   def plot_single(self, ax=None, textsize=20):
      ax = super(idropsPlotFactory, self).plot_single(idrops_markers, ax=ax,
                                                 textsize=textsize)
      return ax

def plot_dropouts_grid(axes_pad=0.1):
   # Plot 4 dropouts on a 2x2 grid of panels
   fig = plt.figure()
   grid = ImageGrid(fig, 111, nrows_ncols=(2, 2), 
                    axes_pad=axes_pad, aspect=False)
   udrops = UdropsPlotFactory()
   udrops.plot_single(ax=grid[0])
   bdrops = BdropsPlotFactory()
   bdrops.plot_single(ax=grid[1])
   vdrops = VdropsPlotFactory()
   vdrops.plot_single(ax=grid[2])
   idrops = idropsPlotFactory()
   idrops.plot_single(ax=grid[3])



def plot_udrops_grid(**grid_kwargs):
   factory = PlotFactory(udrops_catalogs, udrops_field_names, 'f606w', 
                         catalog_dir=udrops_catalog_dir,
                         gfmag_lims=udrops_gfmag_lims)
   grid = factory.plot_grid(**grid_kwargs)
   return grid

# def plot_dropout_single(markers, catalogs=udrops_catalogs, field_names=udrops_field_names, filter_name='f606w', sample='z~3.2', catalog_dir=udrops_catalog_dir, gfmag_lims=udrops_gfmag_lims, colors=udrops_galfit_fc, ax=None, textsize=18):
#    factory = PlotFactory(catalogs, field_names, filter_name, 
#                          catalog_dir=catalog_dir,
#                          gfmag_lims=gfmag_lims)
#    if ax == None:
#       fig = plt.figure()
#       ax = fig.add_subplot(111)
#    ax.set_ylim(0.01, 3.0)
#    ax = factory.plot_single(markers, gfmag_lims, ax=ax, 
#                             colors=colors, sample=sample)
#    ax.legend(loc=3)
#    # Add text for each GALFIT magnitude limit
#    ylims = ax.get_ylim()
#    ## First, UDF-depth
#    ax.text(gfmag_lims[0]-0.1, ylims[1]*0.8, 
#            r'UDF$\rightarrow$', 
#            color=colors[0], size=textsize,
#            ha='right', va='top')
#    ## Then, the shallower depths
#    if len(field_names) > 2:
#       fn = r'%s+%s$\rightarrow$' % (udrops_field_names[1].upper(), udrops_field_names[2].upper())
#    else:
#       fn = field_names[1]
#    ax.text(gfmag_lims[1]-0.1, ylims[1]*0.8, 
#            fn, size=textsize, 
#            ha='right', va='top', color=colors[-1])
#    plt.draw()
#    return ax




   # def collapse_dict(self, input_dict, dtype='float'):
   #    """
   #    Concatenate all numpy arrays stored in input_dict. Assume that input_dict
   #    has keys that are self.field_names.
   #    """
   #    output = np.zeros(0, dtype)
   #    for fn in self.field_names:
   #       output = np.concatenate([output, input_dict[fn]])
   #    return output

   # def collect_attr(self, colname):
   #    x = {}
   #    for fn in self.field_names:
   #       x[fn] = getattr(self.c[fn], colname)
   #    return x

# class PlotMagRe(PlotMeasure):
#    def __init__(self, catalogs, sample_name, field_names, filtername, format='FITS', mag='magout', re='reout_arcsec', mag_err='magout_err', re_err='reout_err_arcsec', goodfit='f850lp_gfflag', yunit='arcsec'):
#       super(PlotMagRe, self).__init__(catalogs, sample_name, field_names, 
#                                       filtername, format='FITS')
#       self.mag = self.collect_attr(mag)
#       self.re = self.collect_attr(re)
#       self.logre = np.where(self.re > 0, np.log10(self.re), -10.)
#       self.logr = self.logre
#       self.mag_err = self.collect_attr(mag_err)
#       self.re_err = self.collect_attr(re_err)
#       self.re_unit = yunit
#       self.plot = self.collect_attr(goodfit)

#    def find_nearest(self, mag0, re0, mag_array, re_array):
#       """
#       Find the nearest point from mag_array and re_array to (mag0, re0)
#       """
#       dist = ((mag_array-mag0)**2 + (re_array-re0)**2)
#       imin = np.argsort(dist)[0]
#       return imin

#    def scatter_mag_logRe(self, mag, Re, marker, size, color, mag_err=[], Re_err=[], ecolor=None, ax=None, ebar_all=False, ebar_one=False, ebar_loc=[25.0,0.3], yunit='arcsec', label="", scatter_kw={}, ebar_kw={}):
#       """
#       Plot GALFIT-measure magnitudes & log(Re) in arcsec.
#       User can decide to add errorbars to all objects (errorbar_all=True), or 
#       add one errorbar to the point closest to a specified (mag, Re) pair (Re
#       in arcsec).
#       However, this is not for plotting the measured values overlaid on top of
#       the RL distribution... The class PlotRLDistFit in plot_fits.py does 
#       that.
#       """
#       if ax == None:
#          fig = plt.figure()
#          ax = fig.add_subplot(111)
#          ax.set_yscale(yscale)
#       if ecolor==None:
#          ecolor = color
#       # start plotting
#       ax.set_yscale('log')
#       ax.scatter(mag, Re, marker=marker, s=size, c=color, label=label, 
#                  **scatter_kw)
#       if ebar_all:
#          assert len(Re_err) > 0
#          ax.errrobar(mag, Re, xerr=mag_err, yerr=Re_err, fmt=None, 
#                      ecolor=ecolor, **ebar_kw)
#       elif ebar_one:
#          assert len(Re_err) > 0
#          mag0, Re0 = ebar_loc
#          imin = self.find_nearest(mag0, Re0, mag, Re)
#          ax.errorbar([mag[imin]], [Re[imin]], xerr=[mag_err[imin]], 
#                      yerr=[Re_err[imin]], fmt=None, ecolor=color, **ebar_kw)

# class PlotLBGMagRe(PlotMagRe):
#    # plots the LBG sample given a fitting parameter file, and creates the 
#    # galaxy_sample.LBGSample instance to initialize the values for mag and Re
#    def __init__(self, paramfile, filtername, sample_name, yunit='arcsec'):
#       LBG = gs.LBGSample(paramfile, filtername)
#       self.sample_name = sample_name
#       self.field_names = LBG.fields
#       self.filtername = filtername
#       self.mag = {}
#       self.re = {}
#       self.logr = {}
#       self.mag_err = {}
#       self.re_err = {}
#       self.plot = {}
#       for f in LBG.fields:
#          sample_f = LBG.galaxy_samples[f]
#          self.mag[f] = sample_f.mag
#          self.logr[f] = sample_f.logr
#          self.re[f] = 10.**sample_f.logr
#          self.plot[f] = np.ones(len(sample_f.mag), 'bool')
#          self.mag_err[f] = sample_f.get_col('mag_err')[sample_f.goodfit]
#          self.re_err[f] = sample_f.get_col('re_err')[sample_f.goodfit]
#       self.logre = self.logr
#       self.titlesize = 20
#       self.axlabelsize = 16
#       self.re_unit = yunit


