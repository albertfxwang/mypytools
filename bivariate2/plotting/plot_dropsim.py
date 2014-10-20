#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from pygoods import Ftable, sextractor
import scipy.stats
from stats import robust


dropsim_dir = '/Users/khuang/Dropbox/Research/bivariate/bivariate_fit/dropsim_catalogs'
udropsim_deep = 'udropsim_goodss_deep_130523.fits'

class PlotSimColors(object):
   def __init__(self, simcatalog):
      self.c = Ftable(simcatalog)

   def get_mag_in(self, b):
      # retrieve input magnitude for a given filter name from the sim catalog
      return getattr(self.c, '%s_mag_in' % b.lower())

   def get_mag_out(self, b, mode):
      # retrieve measured magnitude in a given mode (iso, auto, etc.) for a 
      # given filter from the sim catalog
      return getattr(self.c, '%s_mag_%s' % (b.lower(), mode.lower()))

   def compare_xy(self, x, y, extent, binwd=0.05, **ebar_kwargs):
      # Make a plot with binned statistics laid on top of hexbin density map
      assert len(x) == len(y), "Input x and y should have the same length."
      plt.figure(figsize=(8,6))
      ax = plt.subplot(111)
      ax.hexbin(x, y, extent=extent)
      nbins = int(round((x.max() - x.min()) / binwd))
      xbins = np.linspace(x.min(), x.max(), nbins+1)
      binmed = scipy.stats.binned_statistic(x, y, 'median', bins=xbins)
      binmadn = scipy.stats.binned_statistic(x, y, robust.MADN, bins=xbins)
      ax.errorbar(binmed[1][:-1]+binwd/2, binmed[0], yerr=binmadn[0], fmt='s', 
                  color='yellow', **ebar_kwargs)
      ax.set_xlim(*extent[:2])
      ax.set_ylim(*extent[2:])
      return ax, binmed

   def compare_io_colors(self, band1, band2, colorname, sn_lolim={'H':5.0,'V':5.0,'Y':5.0}, mode='iso', input_extent=None, **ebar_kwargs):
      # Plot output color v.s. input color and calculate the offset between 
      # binned medians
      mag_in_1 = self.get_mag_in(band1)
      mag_in_2 = self.get_mag_in(band2)
      mag_out_1 = self.get_mag_out(band1, mode)
      mag_out_2 = self.get_mag_out(band2, mode)
      color_input = mag_in_1 - mag_in_2
      color_output = mag_out_1 - mag_out_2
      if input_extent == None:
         input_extent = [color_input.min(), color_input.max()]
      extent = list(input_extent) * 2
      # calculate S/N criteria using AUTO magnitudes
      sncrit = (self.c.detect==True)
      for b in sn_lolim.keys():
         flux = getattr(self.c, '%s_flux_auto' % b.lower())
         fluxerr = getattr(self.c, '%s_fluxerr_auto' % b.lower())
         sn = flux / fluxerr
         sncrit = np.logical_and(sncrit, (sn >= sn_lolim[b]))
      sncrit = sncrit & (color_input >= input_extent[0]) & (color_input < input_extent[1])
      print "Number of input simulated sources: %d" % len(self.c.d)
      print "Number of sources in the plot: %d" % np.sum(sncrit)
      color_input = color_input[sncrit]
      color_output = color_output[sncrit]
      ax, binmed = self.compare_xy(color_input, color_output, binwd=0.05, 
                           extent=extent, **ebar_kwargs)
      offset = binmed[0] - (binmed[1][:-1] + binwd/2.)
      ax.text(0.05, 0.9, 'median offset = %.3f' % np.median(offset), 
               ha='left', va='bottom', transform=ax.transAxes,  size='large', 
               bbox=dict(boxstyle='round', facecolor='azure', alpha=0.8))
      ax.set_xlabel('%s (Input)' % colorname.upper())
      ax.set_ylabel('%s (%s)' % (colorname.upper(), mode.upper()))
      ax.set_title(self.c.filename)
      return ax, offset

   def compare_io_mag_colors(self, refband, band1, band2, colorname, sn_lolim={'H':5.0,'V':5.0,'Y':5.0}, mode='iso', input_extent=None, yextent=[-0.5,0.5], relim=[5.,50.], addtext="", **ebar_kwargs):
      # Plot input magnitude in refband v.s. (output - input) color
      mag_in_1 = self.get_mag_in(band1)
      mag_in_2 = self.get_mag_in(band2)
      mag_out_1 = self.get_mag_out(band1, mode)
      mag_out_2 = self.get_mag_out(band2, mode)
      color_input = mag_in_1 - mag_in_2
      color_output = mag_out_1 - mag_out_2
      color_diff = color_output - color_input
      mag_in_ref = self.get_mag_in(refband)
      if input_extent == None:
         input_extent = [mag_in_ref.min(), mag_in_ref.max()]
      extent = input_extent + yextent
      # calculate S/N criteria using AUTO magnitudes
      sncrit = (self.c.detect==True)
      for b in sn_lolim.keys():
         flux = getattr(self.c, '%s_flux_auto' % b.lower())
         fluxerr = getattr(self.c, '%s_fluxerr_auto' % b.lower())
         sn = flux / fluxerr
         sncrit = np.logical_and(sncrit, (sn >= sn_lolim[b]))
      sncrit = sncrit & (mag_in_ref >= input_extent[0]) & (mag_in_ref < input_extent[1])
      # also only select a range of input galaxy size
      recrit = (self.c.re_in>=relim[0]) & (self.c.re_in<relim[1])
      sncrit = sncrit & recrit
      color_diff = color_diff[sncrit]
      mag_in_ref = mag_in_ref[sncrit]
      print "Number of input simulated sources: %d" % len(self.c.d)
      print "Number of sources in the plot: %d" % np.sum(sncrit)
      ax, binmed = self.compare_xy(mag_in_ref, color_diff, binwd=0.2, 
                                   extent=extent, **ebar_kwargs)
      ax.text(0.05, 0.9, 'median offset = %.3f' % np.median(binmed[0]), 
               ha='left', va='bottom', transform=ax.transAxes,  size='large', 
               bbox=dict(boxstyle='round', facecolor='azure', alpha=0.8))
      if addtext:
         ax.text(0.05, 0.05, addtext, ha='left', va='bottom', 
                 transform=ax.transAxes, size='large',
                 bbox=dict(boxstyle='round', facecolor='azure', alpha=0.8))
      ax.set_xlabel('%s magnitude (Input)' % refband.upper())
      ax.set_ylabel('[Output - Input] %s color (%s)' % (colorname.upper(), mode.upper()))
      ax.set_title(self.c.filename)
      return ax, binmed[0]

